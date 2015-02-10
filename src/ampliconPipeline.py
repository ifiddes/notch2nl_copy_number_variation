import sys, os, argparse
import cPickle as pickle

from lib.general_lib import FullPaths, DirType

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import logger, setLoggingFromOptions

import src.models as models


def buildParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq_list", type=str, help="list of fastq files")
    parser.add_argument("--output", "-o", type=DirType, action=FullPaths, default="./amplicon_output/",
            help="base output directory that results will be written to. Default is ./amplicon_output/")
    parser.add_argument("--save_intermediate", action="store_true",
            help="Should we store the intermediates for debugging?")
    return parser


def buildAnalyses(target, output, fastqList, saveInter):
    for fastq in open(fastqList):
        fastq = fastq.rstrip()
        name = os.path.basename(fastq).split("_")[0]
        target.addChildTarget(AmpliconModelWrapper(name, output, fastq, saveInter))


class AmpliconModelWrapper(Target):
    """
    Runs BAM slicer pipeline but without the BAM slicing. Takes local fastq file(s) and runs
    it through all the analyses.
    """
    def __init__(self, uuid, baseOutDir, fastqPath, saveInter):
        Target.__init__(self)
        self.uuid = uuid
        self.baseOutDir = baseOutDir
        self.fastqPath = fastqPath
        self.saveInter = saveInter
        self.outDir = os.path.join(self.baseOutDir, self.uuid)
        #index is a bwa index of the region to be aligned to (one copy of notch2)
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
        if not os.path.exists(bamPath):
            models.alignQuery(self.fastqPath, bamPath, self.getLocalTempDir(), self.uuid, self.index)
        sun = models.FilteredSunModel(self.outDir, self.uuid, bamPath)
        sun.run()
        unfilteredSun = models.UnfilteredSunModel(self.outDir, self.uuid, bamPath)
        unfilteredSun.run()

def main():
    parser = buildParser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    i = Stack(Target.makeTargetFn(buildAnalyses, args=(args.output, args.fastq_list, args.save_intermediate))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == "__main__":
    from src.ampliconPipeline import *
    main()
