"""
Calculates kmer weights for a cohort of individuals.
Expects that there are jellyfish Counts files in the output directory structure.
"""
import sys, os, argparse
import cPickle as pickle
from sonLib.bioio import fastaRead
from lib.general_lib import DirType, FileType
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from collections import Counter
from jobTree.src.bioio import logger, setLoggingFromOptions

#hard coded allele fractions seen in 281 TCGA individuals
avg_c_frac = 1.0 * (2 * 250 + 1 * 18) / 281
avg_d_frac = 1.0 * (2 * 130 + 1 * 121) / 281
avg_frac_dict = {"A":2.0, "B":2.0, "C":avg_c_frac, "D":avg_d_frac, "2":2.0}

def build_parser():
    """
    Builds an argument parser for this run
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", required=True, type=FileType, help="pickled graph to load and add k-mer weights to.")
    parser.add_argument("--new_graph", required=True, type=str, help="where to write the new graph to")
    parser.add_argument("--data_dir", help="output dir full of counts dir/uuid/uuid.Counts.fa", required=True, type=DirType)
    parser.add_argument("--out_dir", help="location to write data to", required=True, type=DirType)
    return parser

class buildDict(Target):
    def __init__(self, uuid, path, out_dir, graph):
        Target.__init__(self)
        self.counts = Counter()
        self.normalizingCounts = Counter()
        self.uuid = uuid
        self.path = path
        self.out_dir = out_dir
        self.graph = graph

    def run(self):
        G = pickle.load(open(self.graph))
        for count, seq in fastaRead(self.path):
            rc = reverseComplement(seq)
            if seq in G.kmers:
                self.counts[seq] += int(count)
            if rc[::-1] != seq and rc in G.kmers:
                self.counts[seq] += int(count)
            if seq in G.normalizingKmers :
                self.normalizingCounts[seq] += int(count)
            if rc[::-1] != seq and rc in G.normalizingKmers:
                self.normalizingCounts[seq] += int(count)
        normalizing = 1.0 * sum(self.normalizingCounts.values()) / len(G.normalizingKmers)
        for kmer in self.counts:
            self.counts[kmer] /= normalizing
        pickle.dump(self.counts, open(os.path.join(self.out_dir, self.uuid + ".counts.pickle"), "w"))        


class Merge(Target):
    def __init__(self, count_files, out_dir, graph, new_graph):
        Target.__init__(self)
        self.count_files = count_files
        self.out_dir = out_dir
        self.graph = graph
        self.new_graph = new_graph

    def run(self):
        dicts = []
        for self.uuid, self.path in self.count_files:
            dicts.append(pickle.load(open(os.path.join(self.out_dir, self.uuid + ".counts.pickle"))))
        counts = reduce(lambda x,y: x+y, dicts)
        with open(os.path.join(self.out_dir, "combined_counts.pickle"), "w") as outf:
            pickle.dump(counts, outf)

        G = pickle.load(open(self.graph))
        weights = {}
        for kmer in G.G.nodes():
            input_sequences = G.G.node[k]['positions'].keys()
            weights[kmer] = (2.0 * counts[kmer])  / (len(self.count_files) * sum(avg_frac_dict[x] for x in input_sequences))
        
        with open(os.path.join(self.out_dir, "weights.pickle"), "w") as outf:
            pickle.dump(weights, outf)
        
        G.weightKmers(weights)
        
        with open(self.new_graph, "w") as outf:
            pickle.dump(G, outf)        


def buildDictWrapper(target, count_files, out_dir, graph, new_graph):
    for uuid, path in count_files:
        target.addChildTarget(buildDict(uuid, path, out_dir, graph))
    target.setFollowOnTarget(Merge(count_files, out_dir, graph, new_graph))


def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    count_files = [[x, os.path.join(args.data_dir, x, x + ".Counts.fa")] for x in os.listdir(args.data_dir)]

    i = Stack(Target.makeTargetFn(buildDictWrapper, args=(count_files, args.out_dir, args.graph, args.new_graph))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from src.find_kmer_weights import *
    main()
