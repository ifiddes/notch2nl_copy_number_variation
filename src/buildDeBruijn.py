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
    parser.add_argument("--genome_counts", "-g", type=argparse.FileType("r"),
        help="Jellyfish kmer count fasta over genome sequences NOT containing region of interest."
        "Counts should be of k-1mers (49mer).")
    parser.add_argument("--offset", "-u", type=int,
        help="Offset from first base of alignment to genome.")
    return parser.parse_args()


def parseJellyfishCounts(file_handle):
    rm = ">\n"
    for count, seq in izip(*[file_handle]*2):
        yield int(count.translate(None, rm)), seq.translate(None, rm)


def main(args):
    args = parse_args(args)
    
    G = DeBruijnGraph(args.kmer_size, args.offset)
    
    for name, seq in fastaRead(args.reference):
        G.addSequences(name, seq)

    for name, seq in fastaRead(args.normalizing):
        G.addNormalizing(name, seq)

    G.finishBuild()
    G.pruneGraph()

    if args.genome_counts is not None:
        G.flagNodes(parseJellyfishCounts(args.genome_counts))
    
    pickle.dump(G, args.out)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
