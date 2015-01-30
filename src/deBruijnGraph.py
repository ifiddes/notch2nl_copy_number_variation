import networkx as nx
from Bio import SeqIO
import logging, string
#I wanted to have this class inherit the nx.DiGraph() directly
#but if I do, then finding subgraphs doesn't work correctly because I
#want to have kmer_size be initialized on instantiation, and then
#the copy function fails. SO has failed to help me so far.
def reverse_complement(s, trans=string.maketrans("ATGC","TACG")):
    return string.translate(s, trans)[::-1]


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
        logging.debug("Initializing a DeBruijnGraph")

        self.kmer_size = kmer_size
        self.G = nx.DiGraph()
        self.has_sequences = False
        self.is_pruned = False
        self.paralogs = set()
        self.kmers = set()
        self.reverse_kmers = set()
        self.normalizing_kmers = set()
        self.reverse_normalizing_kmers = set()

    def paralogs(self):
        return list(self.paralogs)

    def nodes(self):
        return self.G.nodes()

    def node(self):
        return self.G.node


    def add_normalizing(self, seqRecord):
        """
        Adds normalizing kmers to the graph. These kmers are from a region of notch2
        that is after the duplication breakpoint, and so has exactly two copies in everyone.

        """
        k = self.kmer_size - 1

        for i in xrange(len(seqRecord)-k):
            self.normalizing_kmers.add(str(seqRecord.seq[i:i+k]).upper())
            self.reverse_normalizing_kmers.add(reverse_complement(str(seqRecord.seq[i:i+k]).upper()))


    def add_sequences(self, seqRecord):
        """
        Adds k1mers to the graph from a Biopython seqRecord. Edges are built as
        the k1mers are loaded.

        """
        k = self.kmer_size - 1
        name = seqRecord.name
        self.paralogs.add(name)

        for i in xrange(len(seqRecord)-k):
            #left and right k-1mers
            km1L, km1R = str(seqRecord.seq[i:i+k]).upper(), str(seqRecord.seq[i+1:i+k+1]).upper()

            if self.G.has_node(km1L) is not True:
                self.G.add_node(km1L, pos=[i], source=[name], count=1)
            else:
                self.G.node[km1L]["pos"].append(i)
                self.G.node[km1L]["source"].append(name)
                self.G.node[km1L]["count"] += 1
            
            if self.G.has_node(km1R) is not True:
                self.G.add_node(km1R, pos=[], source=[], count=0)

            self.G.add_edge(km1L, km1R)
            self.kmers.add(km1L)
            self.reverse_kmers.add(reverse_complement(km1L))

        #need to count the last kmer also
        self.G.node[km1R]["pos"].append(i+1)
        self.G.node[km1R]["source"].append(name)
        self.G.node[km1R]["count"] += 1
        self.kmers.add(km1R)

        self.has_sequences = True
        logging.debug("{} added to DeBruijnGraph.".format(name))


    def prune_graph(self):
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
        logging.debug("Graph pruned.")


    def weakly_connected_subgraphs(self):
        """
        Yields weakly connected subgraphs and their topolgical sort.

        """
        for subgraph in nx.weakly_connected_component_subgraphs(self.G):
            yield (subgraph, nx.topological_sort(subgraph))

    def flag_nodes(self, kmer_iter):
        """
        Iterates over a kmer_iter and flags nodes as being bad.
        This is used to flag nodes whose kmer is represented elsewhere
        in the genome, so that we won't count it later.
        Therefore, kmer_iter should be an iterable yielding
        sequence strings from a Jellyfish count file.
        Note that the kmer counts should be k-1mers

        """
        for k1mer in kmer_iter:
            if self.G.has_node(k1mer):
                self.G.node[k1mer]['bad'] = True
