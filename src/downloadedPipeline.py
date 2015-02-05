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
from src.models import UnfilteredSunModel, FilteredSunModel
from lib.general_lib import formatRatio, rejectOutliers, FullPaths, DirType

from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from jobTree.src.bioio import system, logger, reverseComplement, setLoggingFromOptions


def buildParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "-o", type=DirType, action=FullPaths, default="./output/",
            help="base output directory that results will be written to. Default is ./output/")
    parser.add_argument("--breakpoint_penalty", type=float, default=23.0,
            help="breakpoint penalty used for ILP model.")
    parser.add_argument("--data_penalty", type=float, default=0.68,
            help="data penalty used for ILP model.")
    parser.add_argument("--graph", type=str, action=FullPaths,
            default="./data/graphs/Notch2NL.pickle")
    return parser


class ModelWrapperLocalFiles(Target):
    """
    Runs BAM slicer pipeline but without the BAM slicing. Takes local fastq file(s) and runs
    it through all the analyses.
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

    def combinedPlot(self, ilpDict, filteredSunDict, unfilteredSunDict, maxPos, offset):
        colors = ["#4D4D4D", "#5DA5DA", "#FAA43A", "#60BD68"]
        #used because the ILP model uses single letter labels
        paraMap = {"Notch2NL-A":"A", "Notch2NL-B":"B", "Notch2NL-C":"C", "Notch2NL-D":"D"}
        xvals = np.array(range(offset, offset + maxPos + 1), dtype="int")
        sortedParalogs = ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D"]
        ax = plt.gca()
        fig, plots = plt.subplots(len(sortedParalogs), sharex=True, sharey=True)
        plt.yticks((0, 1, 2, 3, 4))
        plt.suptitle("kmer-DeBruijn ILP and SUN results")
        for i, (p, para) in enumerate(izip(plots, sortedParalogs)):
            p.axis([offset, offset + maxPos + 1, 0, 4])
            p.fill_between(xvals, ilpDict[para], color=colors[i], alpha=0.7)
            p.set_title("{}".format(para))
            if len(unfilteredSunDict[paraMap[para]]) > 0:
                sunPos, sunVals = zip(*unfilteredSunDict[paraMap[para]])
                p.vlines(np.asarray(sunPos), np.zeros(len(sunPos)), sunVals, color="#F17CB0")
            if len(filteredSunDict[paraMap[para]]) > 0:
                sunPos, sunVals = zip(*filteredSunDict[paraMap[para]])
                p.vlines(np.asarray(sunPos), np.zeros(len(sunPos)), sunVals, color="#763D56")            
        fig.subplots_adjust(hspace=0.5)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False) 
        plt.savefig(os.path.join(self.baseOutDir, self.uuid, self.uuid[:8] + ".png"), format="png")
        plt.close()   

    def run(self):
        sun = FilteredSunModel(self.baseOutDir, self.uuid, self.bamPath)
        sun.run()
        unfilteredSun = UnfilteredSunModel(self.baseOutDir, self.uuid, self.bamPath)
        unfilteredSun.run()
        ilp = IlpModel(self.baseOutDir, self.bpPenalty, self.dataPenalty, self.fastqFile, self.uuid, 
                self.graph, self.getLocalTempDir())
        ilp.run()
        self.combinedPlot(ilp.resultDict, sun.resultDict, unfilteredSun.resultDict, ilp.maxPos, ilp.offset)

def buildAnalyses(target, output, breakpoint_penalty, data_penalty, graph, fastqList):
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