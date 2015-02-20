import sys, os, argparse
import cPickle as pickle

from lib.general_lib import FullPaths, DirType

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import logger, setLoggingFromOptions

import src.models as models


def buildParser():
    parser = argparse.ArgumentParser()
    infiles = parser.add_mutually_exclusive_group()
    infiles.add_argument("--fastq", type=str, help="fastq file")
    infiles.add_argument("--fastq_list", type=str, help="list of fastq files")
    parser.add_argument("--name", type=str, help="name")
    parser.add_argument("--output", "-o", type=DirType, action=FullPaths, default="./output/",
                        help="base output directory that results will be written to. Default is ./output/")
    parser.add_argument("--breakpoint_penalty", type=float, default=20.0,
                        help="breakpoint penalty used for ILP model.")
    parser.add_argument("--data_penalty", type=float, default=0.70,
                        help="data penalty used for ILP model.")
    parser.add_argument("--graph", type=str, action=FullPaths,
                        default="./data/graphs/Notch2NL.pickle")
    parser.add_argument("--save_intermediate", action="store_true",
                        help="Should we store the intermediates for debugging?")
    return parser


def buildAnalyses(target, name, output, breakpoint_penalty, data_penalty, graph, fastqList, saveInter=False):
    for fastq in open(fastqList):
        fastq = fastq.rstrip()
        target.addChildTarget(
            ModelWrapperLocalFile(name, output, breakpoint_penalty, data_penalty, graph, fastq, saveInter))


class ModelWrapperLocalFile(Target):
    """
    Runs BAM slicer pipeline but without the BAM slicing. Takes local fastq file(s) and runs
    it through all the analyses.
    """

    def __init__(self, uuid, baseOutDir, bpPenalty, dataPenalty, graph, fastqPath, saveInter):
        Target.__init__(self)
        self.uuid = uuid[:8]
        self.baseOutDir = baseOutDir
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.graph = graph
        self.fastqPath = fastqPath
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
        else:
            bamPath = os.path.join(self.outDir, self.uuid + ".bam")
        models.alignQuery(self.fastqPath, bamPath, self.getLocalTempDir(), self.uuid, self.index)
        sun = models.FilteredSunModel(self.outDir, self.uuid, bamPath)
        sun.run()
        unfilteredSun = models.UnfilteredSunModel(self.outDir, self.uuid, bamPath)
        unfilteredSun.run()
        ilp = models.IlpModel(self.outDir, self.bpPenalty, self.dataPenalty, self.fastqPath, self.uuid, self.graph,
                              self.getLocalTempDir(), saveInter)
        ilp.run()
        models.combinedPlot(ilp.resultDict, ilp.offsetMap, sun.resultDict, unfilteredSun.resultDict, self.uuid,
                            self.outDir)


def main():
    parser = buildParser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    if args.fastq is not None:
        i = Stack(ModelWrapperLocalFiles(args.name, args.output, args.breakpoint_penalty, args.data_penalty, args.graph,
                                         args.fastq, args.save_intermediate)).startJobTree(args)
    else:
        i = Stack(Target.makeTargetFn(buildAnalyses, args=(
            args.name, args.output, args.breakpoint_penalty, args.data_penalty, args.graph, args.fastq_list,
            args.save_intermediate))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")


if __name__ == "__main__":
    from src.fastqPipeline import *

    main()