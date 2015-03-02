import pulp
from collections import Counter, defaultdict
from itertools import izip
from math import log, exp

from src.abstractIlpSolving import SequenceGraphLpProblem
from jobTree.src.bioio import reverseComplement


class Block(object):
    """
    Represents one block (one weakly connected component)
    Each block may have anywhere from 1 to n paralogs represented in
    it. Initialized by passing in all of the nodes from one DeBruijnGraph
    WCC. Stores a mapping between each paralog this WCC represents
    to the start and stop positions in the source sequence.
    Also stores all of the kmers in this block as a set.

    """

    def __init__(self, subgraph, topo_sorted, min_ploidy=0, max_ploidy=4):
        # a set of all kmers represented by this block
        self.kmers = set()
        self.adjustedCount = None

        for kmer in topo_sorted:
            if 'bad' not in subgraph.node[kmer]:
                self.kmers.add(kmer)

        # one variable for each instance of a input sequence
        self.variables = {}

        #since each node has the same set of sequences, and we have a topological sort, we can 
        #pull down the positions of the first node as the genomic start position
        start_node = subgraph.node[topo_sorted[0]]

        #build variables for each instance
        for para in start_node["positions"]:
            for start in start_node["positions"][para]:
                if len(self.kmers) > 0:
                    self.variables[(para, int(start))] = pulp.LpVariable("{}_{}".format(para, start), 
                                   lowBound=min_ploidy, upBound=max_ploidy, cat="Integer")
                else:
                    self.variables[(para, int(start))] = None


    def __len__(self):
        return len(self.kmers)

    def kmerIter(self):
        """iterator over all kmers in this block"""
        for kmer in self.kmers:
            yield kmer

    def variableIter(self):
        """iterator over variables and related values in this block"""
        for (para, start), variable in self.variables.iteritems():
            yield para, start, variable

    def getVariables(self):
        """returns all LP variables in this block"""
        return self.variables.values()

    def getKmers(self):
        """returns set of all kmers in block"""
        return self.kmers

    def getReverseKmers(self):
        """returns reverse complement kmers in this block"""
        return self.reverseKmers

    def get_count(self):
        """returns counts of data seen in this block"""
        return self.count


class KmerModel(SequenceGraphLpProblem):
    """
    Represents a kmer DeBruijnGraph model to infer copy number of highly
    similar paralogs using ILP.

    Each weakly connected component of the graph represents one 'block'
    which will be assigned one LP variable for each source sequence
    within the WCC (each node all comes from the same sequence(s) by the
    pruning method used).

    normalizing is the average kmer count in a region with known copy number of 2.

    This model can be fine tuned with the following parameters:
    breakpointPenalty: determines how much we allow neighboring variables to differ.
        Larger values restricts breakpoints.
    
    dataPenalty: how much do we trust the data? Larger values locks the variables more
        strongly to the data.

    tightness: How much do we want to favor copy number of defaultPloidy?
    
    """

    def __init__(self, deBruijnGraph, normalizing, breakpointPenalty=15, dataPenalty=1, tightness=1, defaultPloidy=2,
                tightness2=1):
        SequenceGraphLpProblem.__init__(self)
        self.blocks = []
        self.block_map = {x[0]: [] for x in deBruijnGraph.paralogs}
        self.offset_map = {x[0]: int(x[1]) for x in deBruijnGraph.paralogs}

        self.normalizing = normalizing
        self.breakpointPenalty = breakpointPenalty
        self.dataPenalty = dataPenalty
        self.tightness = tightness
        self.tightness2 = tightness2
        self.defaultPloidy = defaultPloidy

        self.buildBlocks(deBruijnGraph)
        self.G = deBruijnGraph

    def buildBlocks(self, deBruijnGraph):
        """
        Builds a ILP kmer model. Input:

        DeBruijnGraph, which is a networkx DeBruijnGraph built over the genome region of interest.

        """
        # make sure the graph has been initialized and pruned
        assert deBruijnGraph.is_pruned and deBruijnGraph.has_sequences

        # build the blocks, don't tie them together yet
        for subgraph, topo_sorted in deBruijnGraph.weaklyConnectedSubgraphs():
            b = Block(subgraph, topo_sorted)
            self.blocks.append(b)

        for block in self.blocks:
            for para, start, variable in block.variableIter():
                self.block_map[para].append([start, variable, block])

        #sort the block maps by start positions
        for para in self.block_map:
            self.block_map[para] = sorted(self.block_map[para], key=lambda x: x[0])

        #now we tie the blocks together
        for para in self.block_map:
            #filter out all blocks without variables (no kmers)
            variables = [v for s, v, b in self.block_map[para] if v is not None]
            for i in xrange(1, len(variables)):
                var_a, var_b = variables[i - 1], variables[i]
                self.constrain_approximately_equal(var_a, var_b, penalty=self.breakpointPenalty)

        #tie each variable to be approximately equal to copy number 2 subject to the tightness constraint
        for block in self.blocks:
            for para, start, variable in block.variableIter():
                self.constrain_approximately_equal(self.defaultPloidy, variable, penalty=self.tightness)

        #now we force all Notch2 variables to be equal to 2
        for s, v, b in self.block_map["Notch2"]:
            self.add_constraint(v == 2)

        #tie the sum of each block to be approximately equal to the number of input sequences subject to tightness2
        #for block in self.blocks:
        #    self.constrain_approximately_equal(sum(block.getVariables()), 2.0 * len(block.getVariables()), 
        #                                        penalty=self.tightness2)

    def introduceData(self, kmerCounts, k1mer_size=49):
        """
        Introduces data to this ILP kmer model. For this, the input is assumed to be a dict 
        representing the results of kmer counting a WGS dataset (format seq:count)

        """
        for block in self.blocks:
            if len(block) > 0:
                count = sum(kmerCounts.get(x, 0) * self.G.G.node[x]['weight'] for x in block.getKmers())
                adjustedCount = (1.0 * count) / (len(block) * self.normalizing)
                block.adjustedCount = adjustedCount

                self.constrain_approximately_equal(adjustedCount, sum(block.getVariables()), penalty=self.dataPenalty)

    def reportCopyMap(self):
        """
        Reports copy number for the solved ILP problem, exploding out so there is a value
        for every x position. Used for graphs.
        """
        copy_map = defaultdict(list)
        for para in self.block_map:
            # find the first variable for this one and extrapolate backwards
            for i in xrange(len(self.block_map[para])):
                start, var, block = self.block_map[para][i]
                if var is not None:
                    prevVal = pulp.value(var)
                    for j in xrange(start):
                        copy_map[para].append(prevVal)
                    break
            for j in xrange(i + 1, len(self.block_map[para])):
                start, var, block = self.block_map[para][j - 1]
                stop, next_var, next_block = self.block_map[para][j]
                if var is None:
                    c = prevVal
                else:
                    c = pulp.value(var)
                    prevVal = c
                for k in xrange(start, stop):
                    copy_map[para].append(c)
            c = pulp.value(next_var)
            for k in xrange(stop, stop + len(next_block)):
                copy_map[para].append(c)
        return copy_map, self.offset_map

    def reportCondensedCopyMap(self):
        """
        Reports copy number for each ILP variable, once.
        format [position, span, value]
        """
        copy_map = defaultdict(list)
        prevVar = self.defaultPloidy
        for para in self.block_map:
            offset = self.offset_map[para]
            for start, var, block in self.block_map[para]:
                if var is not None:
                    copy_map[para].append([start + offset, len(block), pulp.value(var)])
                    prevVar = pulp.value(var)
                else:
                    copy_map[para].append([start + offset, len(block), prevVar])
        return copy_map

    def reportCondensedNormalizedRawDataMap(self):
        """
        Reports the raw counts seen at each variable. This is normalized by the number of variables in the block.
        """
        copy_map = defaultdict(list)
        prevVar = 2
        for para in self.block_map:
            offset = self.offset_map[para]
            for start, var, block in self.block_map[para]:
                if var is not None:
                    copy_map[para].append([start + offset, len(block), block.adjustedCount / len(block.getVariables())])
                    prevVar = block.adjustedCount / len(block.getVariables())
                else:
                    copy_map[para].append([start + offset, len(block), prevVar])
        return copy_map

    def getOffsets(self):
        return self.offset_map