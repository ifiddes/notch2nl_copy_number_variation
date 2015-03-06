import networkx as nx
import string
from collections import Counter, defaultdict
from jobTree.src.bioio import reverseComplement


class DeBruijnGraph(object):
    """
    Represents a DeBruijnGraph using networkx. Sequences are stored on each node as a
    k-1mer with edges pointing towards the next k-1mer in the original sequence.

    When initialized, a kmer size must be provided.

    To add sequences to this graph, pass Biopython SeqRecords to add_sequences().
    Once you have loaded sequences, prune the graph using prune_graph().

    prune_graph() removes all edges that fit into the rules below. This creates linear
    weakly connected components which represent continuous segments of sequence
    that represent one or more paralogs and can be used to infer copy number.

    1) every node with more than 1 outgoing edge has all outgoing edges removed
    2) every node with more than 1 incoming edge has all incoming edges removed
    3) every edge connecting two nodes with a different # of input sequences is removed
    This by definition removes all cycles as well as all strong connected components.

    """

    def __init__(self, kmer_size):
        self.kmer_size = kmer_size
        self.G = nx.Graph()
        self.has_sequences = False
        self.is_pruned = False
        self.paralogs = []
        self.kmers = set()
        self.normalizingKmers = set()
        self.weights = {}

    def nodes(self):
        return self.G.nodes()

    def node(self):
        return self.G.node

    def addNormalizing(self, name, seq):
        """
        Adds normalizing kmers to the graph. These kmers are from a region of notch2
        that is after the duplication breakpoint, and so has exactly two copies in everyone.

        These kmers are stored as whichever strand comes first lexographically.
        This is how jellyfish in the -C mode will report the kmer.

        """
        k = self.kmer_size - 1

        for i in xrange(len(seq) - k):
            s = self.strandless(seq[i:i + k].upper())
            if "N" in s:
                continue
            self.normalizingKmers.add(s)

    def constructNodes(self, name, offset, seq):
        """ 
        constructs left and right nodes for each kmer, adding a sequence edge between nodes.
        """
        self.paralogs.append([name, offset])
        for i in xrange(len(seq) - self.kmer_size + 1):
            kmer = strandless(seq[i:i + self.kmer_size].upper())
            self.kmers.add(kmer)
            if "N" in kmer:
                continue
            l = kmer + "_L"
            r = kmer + "_R"
            # should not be possible to have just left or just right
            if self.G.has_node(l) is not True and self.G.has_node(r) is not True:
                assert not(self.G.has_node(l) or self.G.has_node(r))
                self.G.add_node(l)
                self.G.add_node(r)
                self.G.add_edge(l, r, count=1, positions=defaultdict(list), weight=2.0)
                self.G.edge[l][r]['positions'][name].append(i)
            else:
                self.G.edge[l][r]['count'] += 1
                self.G.edge[l][r]['positions'][name].append(i)
            assert len(self.G) % 2 == 0

    def constructAdjacencies(self, seq):
        """
        Constructs adjacency edges for each the graph.
        """
        prev = seq[:self.kmer_size].upper()
        prev_strandless = strandless(prev)
        for i in xrange(1, len(seq) - self.kmer_size + 1):
            prev_size = len(self.G)
            kmer = seq[i:i + self.kmer_size].upper()
            if "N" in kmer or "N" in prev:
                continue
            kmer_strandless = strandless(kmer)
            if prev == prev_strandless:
                # exiting right side of previous kmer
                if kmer == kmer_strandless:
                    # entering left side of next kmer
                    self.G.add_edge(prev + "_R", kmer + "_L")
                else:
                    # entering right side of next kmer
                    self.G.add_edge(prev + "_R", reverseComplement(kmer) + "_R")
            else:
                # exiting left side of previous kmer
                if kmer == kmer_strandless:
                    # entering left side of next kmer
                    self.G.add_edge(reverseComplement(prev) + "_L", kmer + "_L")
                else:
                    # entering right side of next kmer
                    self.G.add_edge(reverseComplement(prev) + "_L", reverseComplement(kmer) + "_R")
            assert prev_size == len(self.G)
            prev = kmer
            prev_strandless = kmer_strandless

    def pruneGraph(self):
        """
        Creates unitig graphs out of the graph by removing all adjacency edges which fit the following rules:
        1) adjacent nodes have multiple non self-loop edges
        2) adjacent nodes have different counts
        """
        to_remove = []
        for n in self.G.nodes():
            kmer = removeLabel(n)
            adjacency_edges = [(a, b) for a, b in self.G.edges(n) if removeLabel(a) != removeLabel(b)]
            if len(adjacency_edges) > 1:
                # remove all adjacency edges from nodes with more than one adjacency edge
                to_remove.extend(adjacency_edges)
                # remove the adjacency edge if it touches a node with a different count
            else:
                # networkx always reports this node first, right?
                a, b = adjacency_edges[0]
                assert a == n
                if self.G.node[a]['count'] != self.G.node[b]['count']:
                    to_remove.extend(adjacency_edges)
        for a, b in to_remove:
            if self.G.has_edge(a, b) and a != b:
                # don't remove self loop edges
                self.G.remove_edge(a, b)  

    def finishBuild(self, graphviz=False):
        """
        Finishes building the graph.

        If graphviz is true, adds a label tag to each sequence edge to improve understandability in graphviz
        """
        if graphviz is True:
            self_loops = set(self.G.selfloop_edges())
            for edge in self.G.edges():
                if edge not in self_loops and removeLabel(edge[0]) == removeLabel(edge[1]):
                    # sequence edge
                    l = removeLabel(edge[0]) + "\\n" + " - ".join([": ".join([y, ", ".join([str(x) for x in self.G.edge[edge[0]][edge[1]]['positions'][y]])]) for y in self.G.edge[edge[0]][edge[1]]['positions']]) + "\\ncount: " + str(self.G.edge[edge[0]][edge[1]]["count"])
                    self.G.edge[edge[0]][edge[1]]["label"] = l
                else:
                    self.G.edge[edge[0]][edge[1]]['penwidth'] = 2

        self.paralogs = sorted(self.paralogs, key=lambda x: x[0])

        #assert len(self.kmers.intersection(self.normalizingKmers)) == 0

    def connectedComponentIter(self):
        """
        Yields connected components.

        """
        for subgraph in nx.connected_component_subgraphs(self.G):
            yield subgraph

    def flagNodes(self, kmer_iter):
        """
        Iterates over a kmer_iter and flags nodes as being bad.
        This is used to flag nodes whose kmer is represented elsewhere
        in the genome, so that we won't count it later.

        """
        for k in kmer_iter:
            k = k.rstrip()
            assert k in self.kmers
            self.G.edge[k + "_L"][k + "_R"]['bad'] = True

    def weightKmers(self, weightDict):
        """
        Takes a python dictionary mapping k1mers to an empirically derived
        weight. Applies a weight tag to each k1mer in the graph.

        """
        for k, w in weightDict.iteritems():
            assert k in self.kmers
            self.G.edge[k + "_L"][k + "_R"]['weight'] = w

def removeLabel(edge):
    """
    removes the left/right label from a edge, returning the actual kmer
    """
    return edge[:-2]


def strandless(k):
    """
    Returns the strandless version of this kmer. This is defined as whichever comes first, the kmer or the
    reverse complement of the kmer lexicographically.
    """
    return sorted([k, reverseComplement(k)])[0]