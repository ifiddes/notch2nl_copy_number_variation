"""
Calculates kmer weights for a cohort of individuals.
Expects that there are jellyfish Counts files in the output directory structure.
"""
import sys, os, argparse
import cPickle as pickle
from sonLib.bioio import fastaRead
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
    parser.add_argument("--graph", help="pickled deBruijnGraph to be labeled", required=True)
    parser.add_argument("--data_dir", help="output dir full of counts dir/uuid/uuid.Counts.fa", required=True)
    parser.add_argument("--out_dir", help="location to write data to", required=True)
    return parser

def buildDict(target, f, out_dir):
    counts = Counter()
    uuid, path = f
    for count, kmer in fastaRead(path):
        counts[kmer] += int(count)
    pickle.dump(counts, open(os.path.join(out_dir, uuid + ".counts.pickle"), "w"))

def buildDictWrapper(target, count_files, out_dir, graph):
    for f in count_files:
        target.addChildTargetFn(buildDict, args=(f, out_dir))
    target.setFollowOnTargetFn(merge, args=(count_files, out_dir, graph))

def merge(count_files, out_dir):
    dicts = []
    for uuid, path in count_files:
        dicts.append(pickle.load(open(os.path.join(out_dir, uuid + ".counts.pickle"))))
    counts = reduce(lambda x,y: x+y, dicts)
    pickle.dump(final, open(os.path.join(out_dir, "combined_counts.pickle"), "w"))

    weights = {}
    for kmer in G.G.nodes():
        input_sequences = [x.split("_")[0][-1] for x in G.G.node[kmer]['label'].split(", ")]
        population_counts = counts[kmer]
        weight = (1.0 * population_counts) / (len(count_files) * sum(avg_frac_dict[x] for x in input_sequences))
        weights[kmer] = weight
    with open(os.path.join(out_dir, "weights.pickle")) as outf:
        pickle.dump(weights, outf)

    G = pickle.load(open(graph))
    for kmer in G.G.nodes():
        G.G.node[kmer]['weight'] = weights[kmer]
    with open(os.path.join(out_dir, "finished_weighted_graph.pickle"), "w") as outf:
        pickle.dump(G, outf)

def main():
    parser = build_parser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    count_files = [[x, os.path.join("output", x, x + ".Counts.fa")] for x in os.listdir(args.data_dir)]

    i = Stack(Target.makeTargetFn(buildDictWrapper, args=(count_files, args.out_dir, args.graph)))


if __name__ == '__main__':
    from find_kmer_weights import *
    main()
