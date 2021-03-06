import sys, os, argparse
import cPickle as pickle

from lib.general_lib import DirType, FileType

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import logger, setLoggingFromOptions

import src.models as models


def buildParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "-o", type=DirType, default="./tmp_output/",
                        help=("base output directory that results will be written to. Default is ./output/"
                              "This where files will be hunted for."))
    parser.add_argument("--breakpoint_penalty", type=float, default=60.0,
                        help="breakpoint penalty used for ILP model.")
    parser.add_argument("--data_penalty", type=float, default=1.3,
                        help="data penalty used for ILP model.")
    parser.add_argument("--tightness_penalty", type=float, default=0.25,
                        help="How closely should a copy number of 2 be enforced?")
    parser.add_argument("--graph", type=FileType,
                        default="./data/graphs/masked_weighted.pickle")
    parser.add_argument("--kmer_size", type=int, default=49, help="kmer size")
    parser.add_argument("--save_intermediate", action="store_true",
                        help="Should we store the intermediates for debugging?")
    return parser


class ModelWrapperDownloadedFiles(Target):
    """
    Runs the models on all fastq files found in the output folder. Will generate BAMs and counts as necessary.
    """
    def __init__(self, uuid, baseOutDir, bpPenalty, dataPenalty, tightnessPenalty, graph, kmerSize, saveInter):
        Target.__init__(self)
        self.uuid = uuid[:8]
        self.baseOutDir = baseOutDir
        self.outDir = os.path.join(self.baseOutDir, self.uuid)
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.tightness = tightnessPenalty
        self.graph = graph
        self.kmerSize = kmerSize
        self.saveInter = saveInter
        # index is a bwa index of the region to be aligned to (one copy of notch2)
        self.index = "./data/SUN_data/hs_n2.masked.fa"
        if not os.path.exists(self.baseOutDir):
            os.mkdir(self.baseOutDir)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

    def run(self):
        fastqPath = os.path.join(self.outDir, self.uuid + ".fastq")
        if not os.path.exists(fastqPath):
            raise RuntimeError("fastq not in fastqPath. Schema is <uuid>.fastq. uuid: {}".format(self.uuid))
        if self.saveInter is not True:
            bamPath = os.path.join(self.getLocalTempDir(), self.uuid + ".bam")
        else:
            bamPath = os.path.join(self.outDir, self.uuid + ".bam")
        if not os.path.exists(bamPath):
            models.alignQuery(fastqPath, bamPath, self.getLocalTempDir(), self.uuid, self.index)
        sun = models.FilteredSunModel(self.outDir, self.uuid, bamPath)
        sun.run()
        unfilteredSun = models.UnfilteredSunModel(self.outDir, self.uuid, bamPath)
        unfilteredSun.run()
        ilp = models.IlpModel(self.outDir, self.bpPenalty, self.dataPenalty, self.tightness,
                            fastqPath, self.uuid, self.graph, self.getLocalTempDir(), self.kmerSize, self.saveInter, 
                            sun.C, sun.D)
        ilp.run()
        models.combinedPlot(ilp.resultDict, ilp.rawCounts, ilp.offsetMap, unfilteredSun.hg38ResultDict, self.uuid, self.outDir)


def buildAnalyses(target, output, breakpoint_penalty, data_penalty, tightness_penalty, graph, kmerSize, saveInter):
    for uuid in os.listdir(output):
        target.addChildTarget(
            ModelWrapperDownloadedFiles(uuid, output, breakpoint_penalty, data_penalty, tightness_penalty, graph, kmerSize, saveInter))


def main():
    parser = buildParser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(
        args.output, args.breakpoint_penalty, args.data_penalty, args.tightness_penalty, args.graph, args.kmer_size, args.save_intermediate))).startJobTree(
        args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == "__main__":
    from src.downloadedPipeline import *

    main()
