#!/usr/bin/env python2.7
"""

BAM-slicer pipeline.

Takes in a pickled dictionary of cgquery strings and then spins off one jobTree job per genome
in that dictionary. Each genome will undergo SUN analysis on both hg19 and hg38 as well as
Debruijn-Kmer copy number analysis.

"""

import argparse
import cPickle as pickle

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import setLoggingFromOptions, logger

from lib.general_lib import FullPaths, DirType

from src.models import ModelWrapper

def parseArgs():
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

    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)
    return args


def buildAnalyses(target, queries, baseOutDir, bpPenalty, dataPenalty, fullSun, originalSun, Ilp, keyFile):
    for uuid, queryString in queries.iteritems():
        target.addChildTarget(ModelWrapper(uuid, queryString, baseOutDir, bpPenalty, dataPenalty, 
                fullSun, originalSun, Ilp, keyFile))


def main():
    args = parseArgs()
    queries = pickle.load(args.queries)

    if not os.path.exists(args.queries):
        raise RuntimeError("Query file {} does not exist!".format(args.queries))

    if not args.ILP and not args.full and not args.original:
        raise argparse.ArgumentTypeError(("Did not pick any of the three models."
            " Set ILP/full/original flags."))

    if not os.path.exists(args.output):
        logger.info("Making output dir {}".format(args.output))
        os.mkdir(args.output)
    else:
        logger.info("Output dir {} exists".format(args.output))

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(queries, args.output, args.breakpoint_penalty,
            args.data_penalty, args.full, args.original, args.ILP, args.key_file))


if __name__ == "__main__":
    from src.bam_slicer_pipeline import *
    main()