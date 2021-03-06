#!/usr/bin/env python2.7
"""

BAM-slicer pipeline.

Takes in a pickled dictionary of cgquery strings and then spins off one jobTree job per genome
in that dictionary. Each genome will undergo SUN analysis well as Debruijn-Kmer copy number analysis.

"""

import sys, os, argparse
import cPickle as pickle

from lib.general_lib import FullPaths, DirType

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import logger, setLoggingFromOptions

import src.models as models


def buildParser():
    parser = argparse.ArgumentParser()

    parser.add_argument("--queries", "-q", type=argparse.FileType("rb"), default="./queries/queries.pickle",
                        help="pickled query file produced by cgquery_handler.py. Default is ./queries/queries.pickle")
    parser.add_argument("--output", "-o", type=DirType, action=FullPaths, default="./output/",
                        help="base output directory that results will be written to. Default is ./output/")
    parser.add_argument("--breakpoint_penalty", type=float, default=35.0,
                        help="breakpoint penalty used for ILP model.")
    parser.add_argument("--data_penalty", type=float, default=1.0,
                        help="data penalty used for ILP model.")
    parser.add_argument("--tightness_penalty", type=float, default=0.05,
                        help="How closely should a copy number of 2 be enforced?")
    parser.add_argument("--key_file", type=str, action=FullPaths,
                        default="/inside/home/cwilks/haussl_new.key",
                        help="The key file to download protected data from cghub.")
    parser.add_argument("--graph", type=str, action=FullPaths,
                        default="./data/graphs/masked.pickle")
    parser.add_argument("--kmer_size", type=int, default=49, help="kmer size")
    parser.add_argument("--save_intermediate", action="store_true",
                        help="Should we store the intermediates for debugging?")
    return parser


class SlicerModelWrapper(Target):
    """
    This Target runs all of the models.
    First, the fastq is extracted from the BAM slicer via curl and samtools.
    Next, the SUN model is ran (see SunModel)
    Then, the ILP model is ran (see IlpModel)
    Finally, the results of both models is used to build a combined plot.
    """

    def __init__(self, uuid, queryString, baseOutDir, bpPenalty, dataPenalty, tightness, keyFile, graph, kmerSize, saveInter):
        Target.__init__(self)
        self.uuid = uuid[:8]
        self.queryString = queryString
        self.baseOutDir = baseOutDir
        self.outDir = os.path.join(baseOutDir, self.uuid)
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.tightness = tightness
        self.graph = graph
        self.saveInter = saveInter
        self.key = open(keyFile).readline().rstrip()
        self.kmerSize = kmerSize
        self.saveInter = saveInter
        # index is a bwa index of the region to be aligned to (one copy of notch2)
        self.index = "./data/SUN_data/hs_n2.unmasked.fa"
        if not os.path.exists(self.baseOutDir):
            os.mkdir(self.baseOutDir)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

    def run(self):
        if self.saveInter is not True:
            bamPath = os.path.join(self.getLocalTempDir(), self.uuid + ".bam")
            fastqPath = os.path.join(self.getLocalTempDir(), self.uuid + ".fastq")
        else:
            bamPath = os.path.join(self.outDir, self.uuid + ".bam")
            fastqPath = os.path.join(self.outDir, self.uuid + ".fastq")
        models.downloadQuery(fastqPath, self.getLocalTempDir(), self.key, self.queryString, self.uuid)
        models.alignQuery(fastqPath, bamPath, self.getLocalTempDir(), self.uuid, self.index)
        sun = models.FilteredSunModel(self.outDir, self.uuid, bamPath)
        sun.run()
        unfilteredSun = models.UnfilteredSunModel(self.outDir, self.uuid, bamPath)
        unfilteredSun.run()
        ilp = models.IlpModel(self.outDir, self.bpPenalty, self.dataPenalty, self.tightness,
                            fastqPath, self.uuid, self.graph, self.getLocalTempDir(), self.kmerSize, self.saveInter, 
                            sun.C, sun.D)
        ilp.run()
        models.combinedPlot(ilp.resultDict, ilp.rawCounts, ilp.offsetMap, unfilteredSun.hg38ResultDict, self.uuid,
                            self.outDir)


def buildAnalyses(target, queries, baseOutDir, bpPenalty, dataPenalty, tightness, keyFile, graph, kmerSize, saveInter):
    logger.info("Starting to build analyses")
    for uuid, queryString in queries.iteritems():
        target.addChildTarget(SlicerModelWrapper(uuid, queryString, baseOutDir, bpPenalty, dataPenalty, tightness,
                                                 keyFile, graph, kmerSize, saveInter))


def main():
    parser = buildParser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)
    queries = pickle.load(args.queries)

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(queries, args.output, args.breakpoint_penalty,
                                                       args.data_penalty, args.tightness_penalty, args.key_file,
                                                       args.graph, args.kmer_size, args.save_intermediate))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == "__main__":
    from src.bamSlicerPipeline import *

    main()
