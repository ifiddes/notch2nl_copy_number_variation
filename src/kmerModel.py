import pulp
from collections import Counter, defaultdict
import networkx as nx
from itertools import izip, groupby
from operator import itemgetter
from math import log, exp

from src.abstractIlpSolving import SequenceGraphLpProblem
from jobTree.src.bioio import reverseComplement


class Block(object):
    """
    Represents one block (connected component).
    Each block may have anywhere from 1 to n paralogs represented in
    it. Initialized by passing in all of the nodes from one DeBruijnGraph
    CC. Stores a mapping between each paralog this CC represents
    to the start and stop positions in the source sequence.
    Also stores all of the kmers in this block as a set.

    """

    def __init__(self, subgraph, min_ploidy=0, max_ploidy=4):
        self.variables = {}
        self.adjustedCount = None
        # find all sequence edges
        sequence_edges = [x for x in subgraph.edges() if removeLabel(x[0]) == removeLabel(x[1])]
        # find all valid kmers
        self.kmers = {removeLabel(a) for a, b in sequence_edges if 'bad' not in subgraph[a][b]}
        begin_a, begin_b = sequence_edges[0]
        # now we want to find the first and last nodes
        # loop over paralogs
        for p in subgraph[begin_a][begin_b]['positions'].iterkeys():
            # pick a random sequence edge to start at (networkx stores these unordered)
            start_a, start_b = begin_a, begin_b
            # loop over the instances of this paralog
            for i in xrange(len(subgraph[start_a][start_b]['positions'][p])):
                # find the start and stop node for each instance of this paralog
                f = s = subgraph[start_a][start_b]['positions'][p][i]
                # loop over all remaining sequence edges and update the start/stop as necessary
                for new_a, new_b in sequence_edges[1:]:
                    new_s = subgraph[new_a][new_b]['positions'][p][i]
                    new_f = subgraph[new_a][new_b]['positions'][p][i]
                    if new_s < s:
                        s = new_s
                        start_a, start_b = new_a, new_b
                    if new_f > f:
                        f = new_f
                        stop_a, stop_b = new_a, new_b
                assert f - s + 1 >= len(subgraph) / 2
                span = f - s + 1
                if len(self.kmers) > 0:
                    self.variables[(para, s, span)] = pulp.LpVariable("{}_{}".format(para, s), lowBound=min_ploidy, 
                                                                                     upBound=max_ploidy, cat="Integer")
                else:
                    self.variables[(para, s, span)] = None

    def kmerIter(self):
        """iterator over all kmers in this block"""
        for kmer in self.kmers:
            yield kmer

    def variableIter(self):
        """iterator over variables and related values in this block"""
        for (para, start, span), variable in self.variables.iteritems():
            yield para, start, span, variable

    def getVariables(self):
        """returns all LP variables in this block"""
        return [x for x in self.variables.values() if x is not None]

    def getKmers(self):
        """returns set of all kmers in block"""
        return self.kmers


class KmerModel(SequenceGraphLpProblem):
    """
    Represents a kmer DeBruijnGraph model to infer copy number of highly
    similar paralogs using ILP.

    normalizing is the average kmer count in a region with known copy number of 2.

    This model can be fine tuned with the following parameters:
    breakpointPenalty: determines how much we allow neighboring variables to differ.
        Larger values restricts breakpoints.
    
    dataPenalty: how much do we trust the data? Larger values locks the variables more
        strongly to the data.

    tightness: How much do we want to favor copy number of defaultPloidy?
    
    """
    def __init__(self, deBruijnGraph, normalizing, breakpointPenalty=15, dataPenalty=1, tightness=1, inferC=None, 
                 inferD=None, defaultPloidy=2):
        SequenceGraphLpProblem.__init__(self)
        self.blocks = []
        self.block_map = {x[0]: [] for x in deBruijnGraph.paralogs}
        self.offset_map = {x[0]: int(x[1]) for x in deBruijnGraph.paralogs}

        self.normalizing = normalizing
        self.breakpointPenalty = breakpointPenalty
        self.dataPenalty = dataPenalty
        self.tightness = tightness
        self.inferC = inferC
        self.inferD = inferD
        self.defaultPloidy = defaultPloidy

        self.buildBlocks(deBruijnGraph)
        self.G = deBruijnGraph

    def buildBlocks(self, deBruijnGraph):
        """
        Builds a ILP kmer model. Input:

        deBruijnGraph, which is a networkx DeBruijnGraph built over the genome region of interest.

        """
        # build the blocks, don't tie them together yet
        for subgraph in deBruijnGraph.connectedComponentIter():
            b = Block(subgraph)
            self.blocks.append(b)

        for block in self.blocks:
            for para, start, span, variable in block.variableIter():
                self.block_map[para].append([start, span, variable])

        #sort the block maps by start positions
        for para in self.block_map:
            self.block_map[para] = sorted(self.block_map[para], key=lambda x: x[0])

        #now we tie the blocks together
        for para in self.block_map:
            #filter out all blocks without variables (no kmers)
            variables = [var for start, span, var in self.block_map[para] if var is not None]
            for i in xrange(1, len(variables)):
                var_a, var_b = variables[i - 1], variables[i]
                self.constrain_approximately_equal(var_a, var_b, penalty=self.breakpointPenalty)

        #tie each variable to be approximately equal to copy number 2 subject to the tightness constraint
        #for block in self.blocks:
        #    for para, start, variable in block.variableIter():
        #        self.constrain_approximately_equal(self.defaultPloidy, variable, penalty=self.tightness)

        self.exp_dict = {x[0] : self.defaultPloidy for x in deBruijnGraph.paralogs}
        if self.inferC is not None:
            self.exp_dict["Notch2NL-C"] = self.inferC
        if self.inferD is not None:
            self.exp_dict["Notch2NL-D"] = self.inferD
        for block in self.blocks:
            expected = sum(self.exp_dict[p] for p, s, v in block.variableIter() if v is not None)
            self.constrain_approximately_equal(expected, sum(block.getVariables()), penalty=self.tightness)

        #now we force all Notch2 variables to be equal to 2
        if "Notch2" in self.block_map:
            for s, v, b in self.block_map["Notch2"]:
                self.add_constraint(v == 2)

        #if we have previously inferred C/D copy numbers, set those values
        if self.inferC is not None:
            for s, v, b in self.block_map["Notch2NL-C"]:
                self.add_constraint(v == self.inferC)

        if self.inferD is not None:
            for s, v, b in self.block_map["Notch2NL-D"]:
                self.add_constraint(v == self.inferD)

    def introduceData(self, kmerCounts):
        """
        Introduces data to this ILP kmer model. For this, the input is assumed to be a dict 
        representing the results of kmer counting a WGS dataset (format seq:count)

        """
        for block in self.blocks:
            if len(block.kmers) > 0:
                count = sum(kmerCounts[k] * self.G.weights[k] for k in block.getKmers())
                adjustedCount = (1.0 * count) / (len(block.kmers) * self.normalizing)
                block.adjustedCount = adjustedCount
                self.constrain_approximately_equal(adjustedCount, sum(block.getVariables()), penalty=self.dataPenalty)

    def reportCopyMap(self):
        """
        Reports copy number for each ILP variable, once. If a block lacks a variable, reports the previous value.
        format [position, span, value]
        """
        copy_map = defaultdict(list)
        for para in self.block_map:
            prevVar = self.exp_dict[para]
            offset = self.offset_map[para]
            for start, span, var in self.block_map[para]:
                if var is not None:
                    copy_map[para].append([start + offset, span, pulp.value(var)])
                    prevVar = pulp.value(var)
                else:
                    copy_map[para].append([start + offset, span, prevVar])
        return copy_map

    def reportNormalizedRawDataMap(self):
        """
        Reports the raw counts seen at each variable. This is normalized by the number of variables in the block.
        """
        copy_map = defaultdict(list)
        prevVar = 2
        for para in self.block_map:
            offset = self.offset_map[para]
            for start, span, var in self.block_map[para]:
                if var is not None:
                    copy_map[para].append([start + offset, span, block.adjustedCount / len(block.getVariables())])
                    prevVar = block.adjustedCount / len(block.getVariables())
                else:
                    copy_map[para].append([start + offset, span, prevVar])
        return copy_map

    def getOffsets(self):
        return self.offset_map


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