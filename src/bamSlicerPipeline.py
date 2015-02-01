#!/usr/bin/env python2.7
"""

BAM-slicer pipeline.

Takes in a pickled dictionary of cgquery strings and then spins off one jobTree job per genome
in that dictionary. Each genome will undergo SUN analysis on both hg19 and hg38 as well as
Debruijn-Kmer copy number analysis.

"""

import os, argparse
import cPickle as pickle

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, logger

from lib.general_lib import FullPaths, DirType

from src.models import ModelWrapper

def buildParser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--queries", "-q", type=argparse.FileType("rb"), default="./queries/queries.pickle",
            help="pickled query file produced by cgquery_handler.py. Default is ./queries/queries.pickle")
    parser.add_argument("--output", "-o", type=DirType, action=FullPaths, default="./output/",
            help="base output directory that results will be written to. Default is ./output/")
    parser.add_argument("--breakpoint_penalty", type=float, default=15.0,
            help="breakpoint penalty used for ILP model.")
    parser.add_argument("--data_penalty", type=float, default=1.0,
            help="data penalty used for ILP model.")
    parser.add_argument("--no_ILP", action="store_false",
            help="Should the ILP model NOT be ran on these queries?")
    parser.add_argument("--no_full", action="store_false",
            help="Should the SUN model based on the full complement of SUNs NOT be ran on these queries?")
    parser.add_argument("--no_original", action="store_false",
            help="Should the SUN model based on the original SUNs NOT be ran on these queries?")
    parser.add_argument("--key_file", type=str, action=FullPaths,
            default="/inside/home/cwilks/haussl_cghub.key",
            help="The key file to download protected data from cghub.")
    parser.add_argument("--graph", type=str, action=FullPaths,
            default="data/graphs/Notch2NL.pickle")
    return parser


def buildAnalyses(target, queries, baseOutDir, bpPenalty, dataPenalty, fullSun, originalSun, Ilp, 
    keyFile, graph):
    logger.info("Starting to build analyses")
    for uuid, queryString in queries.iteritems():
        target.addChildTarget(ModelWrapper(uuid, queryString, baseOutDir, bpPenalty, dataPenalty, 
                fullSun, originalSun, Ilp, keyFile, graph))


def main():
    parser = buildParser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)
    queries = pickle.load(args.queries)

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(queries, args.output, args.breakpoint_penalty,
            args.data_penalty, args.no_full, args.no_original, args.no_ILP, args.key_file,
            args.graph))).startJobTree(args)


    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == "__main__":
    from src.bamSlicerPipeline import *
    main()