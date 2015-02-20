import sys, os, string, argparse
from itertools import izip
"""
Finds bad kmers. Needs jellyfish counts across an entire genome assembly, and jellyfish across the fasta file
that was used to build the graph.

"""
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome_counts", "-g", required=True, help="Genome counts")
    parser.add_argument("--graph_counts", "-c", required=True, help="Counts of graph region.")
    parser.add_argument("--outfile", "-o", required=True, type=argparse.FileType("w"), help="file to write counts to.")
    return parser.parse_args()


def main():
    args = parse_args()
    rm = ">\n"
    with open(args.graph_counts) as f:
        graph_counts = {x[1].translate(None, rm) : int(x[0].translate(None, rm)) for x in izip(*[f] * 2)}
    print "Starting to go through genome counts"
    bad_count = 0
    good_count = 0
    with open(args.genome_counts) as f:
        for i, (count, seq) in enumerate(izip(*[f] * 2)):
            count = count.translate(None, rm)
            seq = seq.translate(None, rm)
            if seq in graph_counts:
                if graph_counts[seq] < int(count):
                    args.outfile.write(seq + "\n")
                    bad_count += 1
                else:
                    good_count += 1
            if i % 10000000 == 0:
                print "Ran through {} kmers: {} bad, {} good.".format(i, bad_count, good_count)

if __name__ == "__main__":
    main()