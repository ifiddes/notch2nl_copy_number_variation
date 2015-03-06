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
from collections import Counter, defaultdict
from jobTree.src.bioio import logger, setLoggingFromOptions, system
from itertools import izip
import numpy as np

#hard coded allele fractions seen in 201 TCGA individuals
avg_frac_dict = {"Notch2NL-A":2.0, "Notch2NL-B":2.0, "Notch2NL-C":1.843, "Notch2NL-D":0.980, "Notch2":2.0}

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
        self.normalizing = 0
        self.uuid = uuid
        self.path = path
        self.out_dir = out_dir
        self.graph = graph

    def run(self):
        G = pickle.load(open(self.graph))
        normalizing = 0
        with open(self.path) as f:
            rm = ">\n"
            for count, seq in izip(*[f] * 2):
                seq = seq.translate(None, rm)
                count = int(count.translate(None, rm))
                if seq in G.strandlessKmers:
                    self.counts[seq] += count
                elif seq in G.normalizingKmers:
                    self.normalizing += count
        self.normalizing /= len(G.normalizingKmers)
        for kmer in self.counts:
            self.counts[kmer] /= self.normalizing
        pickle.dump(self.counts, open(os.path.join(self.out_dir, self.uuid + ".counts.pickle"), "w"))


class Merge(Target):
    def __init__(self, count_files, out_dir, graph, new_graph):
        Target.__init__(self)
        self.count_files = count_files
        self.out_dir = out_dir
        self.graph = graph
        self.new_graph = new_graph

    def run(self):
        counts = defaultdict(list)
        for d in self.dict_iter():
            for x, y in d.iteritems():
                counts[x].append(y)

        G = pickle.load(open(self.graph))
        kmers = G.kmers

        added_counts = {}
        for k in kmers:
            added_counts[k] = sum(counts[k])

        with open(os.path.join(self.out_dir, "bad_kmers.fasta"), "w") as outf:
            for k in kmers:
                if added_counts[k] == 0:
                    G.G.node[k]['bad'] = True
                    del added_counts[k]
                    outf.write(">{0}\n{0}\n".format(k))

        filtered_kmers = sorted(added_counts.iterkeys())

        with open(os.path.join(self.out_dir, "combined_counts.txt"), "w") as outf:
            for k in filtered_kmers:
                outf.write("{}\t{}\n".format(k, G.weights[k] * added_counts[k]))

        variances = {}
        for k in filtered_kmers:
            variances[k] = np.var(np.asarray(counts[k]))

        with open(os.path.join(self.out_dir, "variances.txt"), "w") as outf:
            for k in filtered_kmers:
                outf.write("{}\t{}\n".format(k, variances[k]))
        
        weights = {}
        for k in filtered_kmers:
            input_sequences = G.G.node[k]['positions'].keys()
            weights[k] = G.weights[k] * (1.0 * len(self.count_files) * sum(avg_frac_dict[x] for x in input_sequences) / (added_counts[k] + 1))

        with open(os.path.join(self.out_dir, "high_weight_bad_kmers.fasta"), "w") as outf:
            for k in filtered_kmers:
                if weights[k] > 4.0:
                    G.G.node[k]['bad'] = True
                    del weights[k]
                    outf.write(">{0}\n{0}\n".format(k))

        with open(os.path.join(self.out_dir, "weights.txt"), "w") as outf:
            for k in weights:
                outf.write("{}\t{}\n".format(k, weights[k]))
        
        G.weightKmers(weights)
        
        with open(self.new_graph, "w") as outf:
            pickle.dump(G, outf)

        system("Rscript src/weights.R {} {} {} {} {}".format(os.path.join(self.out_dir, "combined_counts.txt"), os.path.join(self.out_dir, "weights.txt"), os.path.join(self.out_dir, "variances.txt"), len(self.count_files), "weighting_metrics.pdf"))

    def dict_iter(self):
        for uuid, path in self.count_files:
            yield pickle.load(open(os.path.join(self.out_dir, uuid + ".counts.pickle")))


def strandless(k):
    """
    Returns the strandless version of this kmer. This is defined as whichever comes first, the kmer or the
    reverse complement of the kmer lexicographically.
    """
    return sorted([k, reverseComplement(k)])[0]

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
