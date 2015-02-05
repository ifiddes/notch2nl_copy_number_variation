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
    infiles = parser.add_mutually_exclusive_group()
    infiles.add_argument("--fastq", type=str, help="fastq file")
    infiles.add_argument("--fastq_list", type=argparse.FileType("r"), help="list of fastq files")
    parser.add_argument("--name", type=str, help="name")
    parser.add_argument("--output", "-o", type=DirType, action=FullPaths, default="./output/",
            help="base output directory that results will be written to. Default is ./output/")
    parser.add_argument("--breakpoint_penalty", type=float, default=20.0,
            help="breakpoint penalty used for ILP model.")
    parser.add_argument("--data_penalty", type=float, default=0.70,
            help="data penalty used for ILP model.")
    parser.add_argument("--graph", type=str, action=FullPaths,
            default="./data/graphs/Notch2NL.pickle")
    return parser


class ModelWrapperLocalFiles(Target):
    """
    Runs BAM slicer pipeline but without the BAM slicing. Takes local fastq file(s) and runs
    it through all the analyses.
    """
    def __init__(self, uuid, baseOutDir, bpPenalty, dataPenalty, graph, fastqFile):
        Target.__init__(self)
        self.uuid = uuid
        self.baseOutDir = baseOutDir
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.graph = graph
        #index is a bwa index of the region to be aligned to (one copy of notch2)
        self.index = "./data/SUN_data/hs_n2.masked.fa"
        self.fastqFile = fastqFile

    def alignQuery(self):
        #align the extracted reads to the index
        sortedBamPath = os.path.join(self.getLocalTempDir(), "{}.sorted".format(self.uuid))
        system("bwa mem -v 1 {} {} | samtools view -F 4  -bS - | samtools sort - {}".format(self.index, self.fastqFile, sortedBamPath))
        #samtools appends .bam to sorted bam files
        sortedBamPath += ".bam"
        #filter the SAM records and find the site coverage at each locus, creating VCFs
        remappedBamPath = os.path.join(self.baseOutDir, self.uuid, "{}.remapped.sorted.bam".format(self.uuid))
        header = {"HD": {"VN": "1.3"}, "SQ": [{"LN": 248956422, "SN": "chr1"}]}
        outfile = pysam.Samfile(remappedBamPath, "wb", header=header)
        bamfile = pysam.Samfile(sortedBamPath, "rb")  
        for record in bamfile:
            chrom, span = bamfile.getrname(record.tid).split(":")
            start, end = map(int, span.split("-"))
            record.pos = record.pos + start - 1
            outfile.write(record)
        outfile.close()
        system("samtools index {}".format(remappedBamPath))
        smallerFastqFile = os.path.join(self.baseOutDir, self.uuid, self.uuid[:8] + ".fastq")
        system("samtools view {} | bamstools bamshuf -Ou /dev/stdin {} | samtools bam2fq /dev/stdin > {}".format(remappedBamPath, os.path.join(self.getLocalTempDir(), "tmp"), smallerFastqFile))
        return remappedBamPath, self.fastqFile

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
        if not os.path.exists(os.path.join(self.baseOutDir, self.uuid)):
            os.mkdir(os.path.join(self.baseOutDir, self.uuid))
        bamPath, fastqFile = self.alignQuery()
        sun = FilteredSunModel(self.baseOutDir, self.uuid, bamPath)
        sun.run()
        unfilteredSun = UnfilteredSunModel(self.baseOutDir, self.uuid, bamPath)
        unfilteredSun.run()
        ilp = IlpModel(self.baseOutDir, self.bpPenalty, self.dataPenalty, fastqFile, self.uuid, 
                self.graph, self.getLocalTempDir())
        ilp.run()
        self.combinedPlot(ilp.resultDict, sun.resultDict, unfilteredSun.resultDict, ilp.maxPos, ilp.offset)

def buildAnalyses(target, name, output, breakpoint_penalty, data_penalty, graph, fastqList):
    for fastq in fastqList:
        fastq = fastq.rstrip()
        target.addChildTarget(ModelWrapperLocalFiles(name, output, breakpoint_penalty, data_penalty, graph, fastq))

def main():
    parser = buildParser()
    Stack.addJobTreeOptions(parser)
    args = parser.parse_args()
    setLoggingFromOptions(args)

    if os.path.exists(args.jobTree):
        shutil.rmtree(args.jobTree)

    if args.fastq is not None:
        i = Stack(ModelWrapperLocalFiles(args.name, args.output, args.breakpoint_penalty, args.data_penalty, args.graph, args.fastq)).startJobTree(args)
    else:
        i = Stack(Target.makeTargetFn(buildAnalyses, args=(args.name, args.output, args.breakpoint_penalty, args.data_penalty, args.graph, args.fastq_list))).startJobTree(args)

    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == "__main__":
    from src.fastqPipeline import *
    main()