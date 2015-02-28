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
        self.G = nx.DiGraph()
        self.has_sequences = False
        self.is_pruned = False
        self.paralogs = []
        self.kmers = set()
        self.normalizingKmers = set()

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

        for i in xrange(len(seq) - k):
            s = seq[i:i + k].upper()
            if "N" in s:
                continue
            self.normalizingKmers.add(s)


    def addSequences(self, name, offset, seq):
        """

        Adds k1mers to the graph. Edges are built as the k1mers are loaded.

        """
        k = self.kmer_size - 1
        self.paralogs.append([name, offset])
        paralogPosition = Counter()
        prev = ""
        for i in xrange(len(seq) - k + 1):
            k1mer = seq[i:i + k].upper()
            if "N" in prev or "N" in k1mer:
                paralogPosition[name] += 1
                continue
            else:
                if self.G.has_node(k1mer) is not True:
                    self.G.add_node(k1mer, count=1, positions=defaultdict(list))
                    self.G.node[k1mer]['positions'][name].append(paralogPosition[name])
                    paralogPosition[name] += 1
                else:
                    self.G.node[k1mer]['positions'][name].append(paralogPosition[name])
                    self.G.node[k1mer]['count'] += 1
                    paralogPosition[name] += 1
                if prev is not "":
                    self.G.add_edge(prev, k1mer)
                self.kmers.add(k1mer)
            prev = k1mer

        self.has_sequences = True

    def finishBuild(self, graphviz=False):
        if graphviz is True:
            for node in self.G.nodes():
                l = node + "\\n" + " - ".join([": ".join([y, ", ".join([str(x) for x in self.G.node[node]['positions'][y]])]) for y in self.G.node[node]['positions']]) + "\\ncount: " + str(self.G.node[node]["count"])
                self.G.node[node]["label"] = l

        self.paralogs = sorted(self.paralogs, key=lambda x: x[0])
        for n in self.G.nodes():
                self.G.node[n]['weight'] = 2

        assert len(self.kmers.intersection(self.normalizingKmers)) == 0

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

    def flagNodes(self, kmer_iter):
        """
        Iterates over a kmer_iter and flags nodes as being bad.
        This is used to flag nodes whose kmer is represented elsewhere
        in the genome, so that we won't count it later.
        Note that the kmer counts should be k-1mers

        """
        for k1mer in kmer_iter:
            k1mer = k1mer.rstrip()
            if k1mer in self.kmers:
                self.G.node[k1mer]['bad'] = True

    def weightKmers(self, weightDict):
        """
        Takes a python dictionary mapping k1mers to an empirically derived
        weight. Applies a weight tag to each k1mer in the graph.
        """
        for k1mer, weight in weightDict.iteritems():
            if k1mer in self.kmers:
                self.G.node[k1mer]['weight'] = weight