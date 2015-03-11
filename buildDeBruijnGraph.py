#!/usr/bin/env python2.7
import argparse, sys, os
from itertools import izip
import cPickle as pickle

from jobTree.src.bioio import fastaRead
from src.deBruijnGraph import DeBruijnGraph

"""
Builds a kmer DeBruijnGraph with nodes flagged for kmers present in genome.
Serializes this to disk.

"""

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", "-r", type=str, required=True, 
        help="Reference fasta file")
    parser.add_argument("--normalizing", "-n", type=str, required=True,
        help="Normalizing region fasta file.")
    parser.add_argument("--out", "-o", type=argparse.FileType("wb"), required=True,
        help="File to write pickled DeBruijnGraph to.")
    parser.add_argument("--kmer_size", "-k", type=int, default=50, 
        help="kmer size. Default=50")
    parser.add_argument("--bad_kmers", "-b", type=argparse.FileType("r"),
        help="Text file of kmers that are to be flagged as bad. See findBadKmers.py")
    parser.add_argument("--weights", "-w", type=str,
        help="pickled python dictionary representing empirically derived per-kmer weights.")
    return parser.parse_args()


def main(args):
    args = parse_args(args)
    
    G = DeBruijnGraph(args.kmer_size)
    
    # first pass adds nodes
    for name, seq in fastaRead(args.reference):
        name, offset = name.split("_")[:2]
        G.constructNodes(name, offset, seq)
    # second pass constructs adjacenices
    for name, seq in fastaRead(args.reference):
        G.constructAdjacencies(seq)

    for name, seq in fastaRead(args.normalizing):
        G.addNormalizing(name, seq)
        
    if args.bad_kmers is not None:
        G.flagNodes(args.bad_kmers)

    if args.weights is not None:
        with open(args.weights) as f:
            G.weightKmers(pickle.load(f))

    G.finishBuild()
    G.pruneGraph()
    
    pickle.dump(G, args.out)


if __name__ == '__main__':
    sys.exit(main(sys.argv))

