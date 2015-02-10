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
        #a set of all kmers represented by this block
        self.kmers = set()
        self.reverseKmers = set()

        for kmer in topo_sorted:
            if 'bad' not in subgraph.node[kmer]:
                self.kmers.add(kmer)
                self.reverseKmers.add(reverseComplement(kmer))

        #one variable for each instance of a input sequence
        self.variables = {}

        #since each node has the same set of sequences, and we have a topological sort, we can 
        #pull down the positions of the first node only
        start_node = subgraph.node[topo_sorted[0]]
        
        #build variables for each instance
        for para_start in start_node["label"].split(", "):
            para, start = para_start.split("_")
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
    
    sizeAdjust: adjusts dataPenalty so that more weight is given to larger blocks.
        This reduces noise due to sequencing depth. This is plugged into a exponential
        fit so that very large blocks are not linearly more weighted than medium blocks.
    
    singleVariableWeight: adjusts dataPenalty so that more weight is given to blocks
        that contain only one variable. This block contains SUN positions.
    
    tightness: How much do we want to favor copy number of 2?
    
    """
    def __init__(self, deBruijnGraph, normalizing, breakpointPenalty=15, dataPenalty=1, 
                sizeAdjust=0.20, singleVariableWeight=5):
        SequenceGraphLpProblem.__init__(self)
        self.blocks = []
        self.block_map = { x[0] : [] for x in deBruijnGraph.paralogs }
        self.offset_map = { x[0] : int(x[1]) for x in deBruijnGraph.paralogs }

        self.normalizing = normalizing
        self.breakpointPenalty = breakpointPenalty
        self.dataPenalty = dataPenalty
        self.sizeAdjust = sizeAdjust
        self.singleVariableWeight = singleVariableWeight

        self.buildBlocks(deBruijnGraph)

    def buildBlocks(self, deBruijnGraph):
        """
        Builds a ILP kmer model. Input:

        DeBruijnGraph, which is a networkx DeBruijnGraph built over the genome region of interest.

        """
        #make sure the graph has been initialized and pruned
        assert deBruijnGraph.is_pruned and deBruijnGraph.has_sequences

        #build the blocks, don't tie them together yet
        for subgraph, topo_sorted in deBruijnGraph.weaklyConnectedSubgraphs():
            b = Block(subgraph, topo_sorted)
            self.blocks.append(b)

        for block in self.blocks:
            for para, start, variable in block.variableIter():
                self.block_map[para].append([start, variable, block])

        #sort the block maps by start positions
        for para in self.block_map:
            self.block_map[para] = sorted(self.block_map[para], key = lambda x: x[0])
        
        #now we tie the blocks together
        for para in self.block_map:
            #filter out all blocks without variables (no kmers)
            variables = [v for s, v, b in self.block_map[para] if v is not None]
            for i in xrange(1, len(variables)):
                var_a, var_b = variables[i-1], variables[i]
                self.constrain_approximately_equal(var_a, var_b, penalty=self.breakpointPenalty)

    def introduceData(self, kmerCounts, k1mer_size=49):
        """
        Introduces data to this ILP kmer model. For this, the input is assumed to be a dict 
        representing the results of kmer counting a WGS dataset (format seq:count)

        """

        for block in self.blocks:
            if len(block) > 0:
                count = sum( kmerCounts.get(x, 0) for x in block.getKmers() )
                count += sum( kmerCounts.get(x, 0) for x in block.getReverseKmers() )

                adjustedCount = 2.0 * count / ( len(block) * self.normalizing )

                #weight larger blocks more (data less noisy)
                dp = self.dataPenalty * log(exp(1) + len(block) * self.sizeAdjust)

                #give larger weight to blocks with only one variable
                #since this block contains unique sequence
                if len(block.getVariables()) == 1:
                    dp *= self.singleVariableWeight

                self.constrain_approximately_equal(adjustedCount, sum(block.getVariables()), 
                        penalty=dp)
                
    def reportCopyMap(self):
        """
        Reports copy number for the solved ILP problem, exploding out so there is a value
        for every x position. Used for graphs.
        """
        copy_map = defaultdict(list)
        for para in self.block_map:
            #find the first variable for this one and extrapolate backwards
            for i in xrange(len(self.block_map[para])):
                start, var, block = self.block_map[para][i]
                if var is not None:
                    prevVal = pulp.value(var)
                    for j in xrange(start):
                        copy_map[para].append(prevVal)
                    break
            for j in xrange(i+1, len(self.block_map[para])):
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
        for para in self.block_map:
            offset = self.offset_map[para]
            for start, var, block in self.block_map[para]:
                copy_map[para].append([start + offset, len(block), pulp.value(var)])
        return copy_map