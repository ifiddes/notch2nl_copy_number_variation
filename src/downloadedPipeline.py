import sys, os, pysam, vcf, string, gzip, argparse, shutil
import cPickle as pickle
from itertools import izip
from collections import defaultdict, Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pylab import setp

from src.kmerModel import KmerModel
from src.models import UnfilteredSunModel, FilteredSunModel, IlpModel
from lib.general_lib import formatRatio, rejectOutliers, FullPaths, DirType

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import system, logger, reverseComplement, setLoggingFromOptions


def buildParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "-o", type=DirType, action=FullPaths, default="./output/",
            help=("base output directory that results will be written to. Default is ./output/"
            "For this model is where files will be hunted for."))
    parser.add_argument("--breakpoint_penalty", type=float, default=23.0,
            help="breakpoint penalty used for ILP model.")
    parser.add_argument("--data_penalty", type=float, default=0.68,
            help="data penalty used for ILP model.")
    parser.add_argument("--graph", type=str, action=FullPaths,
            default="./data/graphs/Notch2NL.pickle")
    parser.add_argument("--save_intermediate", action="store_true",
            help="Should we store the intermediates for debugging?")
    return parser


class ModelWrapperDownloadedFiles(Target):
    """
    Runs the models on all fastq files found in the output folder. Will generate BAMs and counts as necessary.
    """
    def __init__(self, uuid, baseOutDir, bpPenalty, dataPenalty, graph):
        Target.__init__(self)
        self.uuid = uuid[:8]
        self.baseOutDir = baseOutDir
        self.outDir = os.path.join(self.baseOutDir, self.uuid)
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.graph = graph
        #index is a bwa index of the region to be aligned to (one copy of notch2)
        self.index = "./data/SUN_data/hs_n2.masked.fa"
        self.fastqFile = os.path.join(self.outDir, self.uuid + ".fastq")
        self.bamPath = os.path.join(self.outDir, self.uuid + ".remapped.sorted.bam")
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
        if not os.path.exists(bamPath):
            models.alignQuery(fastqPath, bamPath, self.getLocalTempDir(), self.queryString, self.uuid, self.index)
        sun = models.FilteredSunModel(self.outDir, self.uuid, bamPath)
        sun.run()
        unfilteredSun = UnfilteredSunModel(self.outDir, self.uuid, bamPath)
        unfilteredSun.run()
        ilp = models.IlpModel(self.outDir, self.bpPenalty, self.dataPenalty, fastqFile, self.uuid, self.graph, self.getLocalTempDir(), saveInter)
        ilp.run()
        models.combinedPlot(ilp.resultDict, sun.resultDict, unfilteredSun.resultDict, ilp.maxPos, ilp.offset, self.uuid, self.outDir)


def buildAnalyses(target, output, breakpoint_penalty, data_penalty, graph):
    for uuid in os.listdir(output):
        target.addChildTarget(ModelWrapperLocalFiles(uuid, output, breakpoint_penalty, data_penalty, graph))


def main():
    parser = buildParser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(args.output, args.breakpoint_penalty, args.data_penalty, args.graph))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == "__main__":
    from src.downloadedPipeline import *
    main()