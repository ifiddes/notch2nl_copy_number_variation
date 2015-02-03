import networkx as nx
import string
from collections import Counter
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
    def __init__(self, kmer_size, offset):

        self.offset = offset
        self.kmer_size = kmer_size
        self.G = nx.DiGraph()
        self.has_sequences = False
        self.is_pruned = False
        self.paralogs = set()
        self.kmers = set()
        self.reverseKmers = set()
        self.normalizingKmers = set()
        self.reverseNormalizingKmers = set()

    def paralogs(self):
        return list(self.paralogs)

    def nodes(self):
        return self.G.nodes()

    def node(self):
        return self.G.node


    def addNormalizing(self, name, seq):
        """
        Adds normalizing kmers to the graph. These kmers are from a region of notch2
        that is after the duplication breakpoint, and so has exactly two copies in everyone.

        """
        k = self.kmer_size - 1

        for i in xrange(len(seq)-k):
            self.normalizingKmers.add(seq[i:i+k].upper())
            self.reverseNormalizingKmers.add(reverseComplement(seq[i:i+k].upper()))


    def addSequences(self, name, seq):
        """

        Adds k1mers to the graph.Edges are built as the k1mers are loaded.

        """
        k = self.kmer_size - 1
        self.paralogs.add(name)
        #TODO - handle names longer than 1 character by truncating
        paralogNodeCount = Counter()

        for i in xrange(len(seq)-k):
            #left and right k-1mers 
            km1L, km1R = seq[i:i+k].upper(), seq[i+1:i+k+1].upper()

            if self.G.has_node(km1L) is not True:
                self.G.add_node(km1L, label=["{}_{}".format(name, paralogNodeCount[name])], count=1)
                paralogNodeCount[name] += 1
            else:
                self.G.node[km1L]["label"].append("{}_{}".format(name, paralogNodeCount[name]))
                self.G.node[km1L]["count"] += 1
                paralogNodeCount[name] += 1
            
            if self.G.has_node(km1R) is not True:
                self.G.add_node(km1R, label=[], count=0)

            self.G.add_edge(km1L, km1R)
            self.kmers.add(km1L)
            self.reverseKmers.add(reverseComplement(km1L))

        #need to count the last kmer also
        self.G.node[km1R]["label"].append("{}_{}".format(name, paralogNodeCount[name]))
        self.G.node[km1R]["count"] += 1
        self.kmers.add(km1R)
        self.reverseKmers.add(reverseComplement(km1R))

        self.has_sequences = True

    def finishBuild(self):
        for node in self.G.nodes():
            self.G.node[node]["label"] = ", ".join(sorted(self.G.node[node]["label"]))

    def pruneGraph(self):
        """
        For each node, if has more than one outgoing edge remove all outgoing edges.
        Do the same for incoming edges. Also, remove all edges between nodes with
        different counts.

        """
        for n in self.G.nodes_iter():
            if len(self.G.in_edges(n)) > 1:
                for n1, n2 in self.G.in_edges(n):
                    self.G.remove_edge(n1, n2)

            if len(self.G.out_edges(n)) > 1:
                for n1, n2 in self.G.out_edges(n):
                    self.G.remove_edge(n1, n2)

            for n1, n2 in self.G.out_edges(n):
                if self.G.node[n1]["count"] != self.G.node[n2]["count"]:
                    self.G.remove_edge(n1, n2)

        self.is_pruned = True


    def weaklyConnectedSubgraphs(self):
        """
        Yields weakly connected subgraphs and their topolgical sort.

        """
        for subgraph in nx.weakly_connected_component_subgraphs(self.G):
            yield (subgraph, nx.topological_sort(subgraph))

    def flagNodes(self, kmer_iter, cutoff=1):
        """
        Iterates over a kmer_iter and flags nodes as being bad.
        This is used to flag nodes whose kmer is represented elsewhere
        in the genome, so that we won't count it later.
        Therefore, kmer_iter should be an iterable yielding
        sequence strings from a Jellyfish count file.
        Note that the kmer counts should be k-1mers

        Cutoff determines the minimum number of counts seen elsewhere
        in the genome to count as a flagged node.

        """
        for count, k1mer in kmer_iter:
            if count >= cutoff:
                if k1mer in self.kmers:
                    self.G.node[k1mer]['bad'] = True
                elif reverseComplement(k1mer) in self.reverseKmers:
                    self.G.node[k1mer]['bad'] = True
