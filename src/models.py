import sys, os, pysam, vcf, string, gzip, math
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
from src.sunIlpModel import SunIlpModel
from lib.general_lib import formatRatio

from jobTree.scriptTree.target import Target
from jobTree.src.bioio import system, logger, reverseComplement, fastaRead


class SunModel(object):
    """
    Runs the SUN model, where alelle fraction at unique sites is used to infer copy number
    """

    def findSiteCoverages(self, bamIn, minDepth=50):
        """
        Runs mpileup on each SUN position and finds allele fractions
        """
        bases = set(["A", "T", "G", "C", "a", "t", "g", "c"])
        resultDict = defaultdict(list)
        nVals = []
        for pos, (para, ref, alt, hg38_pos) in self.wl.iteritems():
            posStr = "chr1:{0}-{0}".format(pos)
            pileUp = pysam.mpileup("-q", "20", "-r", posStr, bamIn)
            if len(pileUp) == 0:
                continue
            pileUpStr = pileUp[0].split()
            if len(pileUpStr) != 6:
                continue
            pileUpResult = Counter(x.upper() for x in pileUpStr[4] if x in bases)
            if ref not in pileUpResult or alt not in pileUpResult:
                continue
            if sum(pileUpResult.itervalues()) < minDepth:
                continue
            frac = formatRatio(pileUpResult[alt], sum(pileUpResult.values()))
            # invert fraction for Notch2 paralogs
            if para == "N":
                frac = 1 - frac
                nVals.append(frac)
            else:
                resultDict[para].append([pos, frac])
        norm = self.findNormalizingFactor(nVals)
        for para, result in resultDict.iteritems():
            resultDict[para] = [[x, 2.0 * y / norm] for x, y in result]
        return resultDict

    def findNormalizingFactor(self, r):
        """
        The values for Notch2NL-A-D are normalized by Notch2. Notch2 does not have CNV.
        """
        return np.mean(np.asarray(r))

    def plotHistograms(self, resultDict):
        path = os.path.join(self.outDir, "{}.{}.png".format(self.uuid, self.__class__.__name__))
        f, plts = plt.subplots(4, sharex=True)
        for para, p in izip(sorted(resultDict.keys()), plts):
            p.set_title("Notch2NL-{}".format(para))
            p.hist(np.asarray(zip(*resultDict[para])[1]), bins=40, range=(0.0, 4.0), color='#1e90ff', normed=True)
            p.yaxis.tick_right()
            setp(p.get_yticklabels(), fontsize=9)
        plt.xticks((0, 1, 2, 3, 4))
        f.subplots_adjust(hspace=0.4)
        plt.savefig(path)

    def convertResultDict(self):
        converted = defaultdict(list)
        for para in self.resultDict:
            for pos, frac in self.resultDict[para]:
                converted[para].append([self.wl[pos][3], frac])
        return converted

    def makeHg19Bedgraphs(self, resultDict):
        if not os.path.exists(os.path.join(self.outDir, "tracks")):
            os.mkdir(os.path.join(self.outDir, "tracks"))
        for para in resultDict:
            path = os.path.join(self.outDir, "tracks", "{}.Notch2NL-{}.{}.hg19.bedGraph".format(
                self.uuid, para, self.__class__.__name__))
            bedHeader = ("track type=bedGraph name={} autoScale=off visibility=full alwaysZero=on yLineMark=2 "
                        "viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n")
            with open(path, "w") as outf:
                outf.write(bedHeader.format(self.uuid + "_" + para))
                for pos, frac in sorted(resultDict[para], key = lambda x:x[0]):
                    outf.write("\t".join(map(str, ["chr1", pos, pos + 1, frac])) + "\n")

    def makeHg38Bedgraphs(self, resultDict):
        if not os.path.exists(os.path.join(self.outDir, "tracks")):
            os.mkdir(os.path.join(self.outDir, "tracks"))
        path = os.path.join(self.outDir, "tracks", "{}.{}.hg38.bedGraph".format(self.uuid, self.__class__.__name__))
        tmp = sorted([y for x in resultDict.itervalues() for y in x], key=lambda x: x[0])
        with open(path, "w") as outf:
            bedHeader = ("track type=bedGraph name={} autoScale=off visibility=full alwaysZero=on yLineMark=2 "
                        "viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n")
            outf.write(bedHeader.format(self.uuid))
            for pos, frac in tmp:
                outf.write("\t".join(map(str, ["chr1", pos, pos + 1, frac])) + "\n")

    def callIlp(self, windowSize=5000, stepSize=250, dataPenalty=1.50, breakpointPenalty=0.70,
                deletionPenalty=2):
        Avals = sorted(self.hg38ResultDict["A"], key = lambda x: x[0])
        Bvals = sorted(self.hg38ResultDict["B"], key = lambda x: x[0])
        model = SunIlpModel(Avals, Bvals, windowSize, stepSize, breakpointPenalty, dataPenalty, deletionPenalty)
        results = model.run()
        if not os.path.exists(os.path.join(self.outDir, "tracks")):
            os.mkdir(os.path.join(self.outDir, "tracks"))
        path = os.path.join(self.outDir, "tracks",
                            "{}.{}.hg38.SUN_ILP.wiggle".format(self.uuid, self.__class__.__name__))
        with open(path, "w") as outf:
            bedHeader = ("track type=wiggle_0 name={} autoScale=off visibility=full alwaysZero=on yLineMark=2 "
                        "viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n")
            outf.write(bedHeader.format(self.uuid))
            for pos, val in results:
                mid = pos + windowSize / 2
                outf.write("variableStep chrom=chr1 span={}\n{} {}\n".format(stepSize, mid, val))

    def inferCopyNumber(self):
        """
        Infer copy number of Notch2NL-C and Notch2NL-D

        """
        C = [x[1] for x in self.resultDict["C"]]
        D = [x[1] for x in self.resultDict["D"]]
        return int(round(sum(C) / len(C), 0)), int(round(sum(D) / len(D), 0))

    def run(self):
        self.resultDict = self.findSiteCoverages(self.bamPath)
        self.C, self.D = self.inferCopyNumber()
        # self.plotHistograms(self.resultDict)
        self.makeHg19Bedgraphs(self.resultDict)
        self.hg38ResultDict = self.convertResultDict()
        self.makeHg38Bedgraphs(self.hg38ResultDict)
        if self.doIlp is True:
            self.callIlp()


class UnfilteredSunModel(SunModel):
    def __init__(self, outDir, uuid, bamPath, doIlp=False):
        self.uuid = uuid
        self.bamPath = bamPath
        self.outDir = outDir
        self.doIlp = doIlp
        # whitelist is a text file of whitelisted SUN positions - in this case, unfiltered
        with open("./data/SUN_data/hg38_unfiltered_whitelist.txt") as wl:
            wl_list = [x.split() for x in wl if not x.startswith("#")]
        # dict mapping genome positions to which paralog has a SUN at that position
        # [paralog, hg19_pos, ref, alt, hg38_pos]
        self.wl = {int(x[1]): [x[0], x[2], x[3], int(x[4])] for x in wl_list}

    def run(self):
        SunModel.run(self)


class FilteredSunModel(SunModel):
    def __init__(self, outDir, uuid, bamPath, doIlp=False):
        self.uuid = uuid
        self.bamPath = bamPath
        self.outDir = outDir
        self.doIlp = doIlp
        # whitelist is a text file of whitelisted SUN positions - in this case, unfiltered
        with open("./data/SUN_data/hg38_whitelist.txt") as wl:
            wl_list = [x.split() for x in wl if not x.startswith("#")]
        # dict mapping genome positions to which paralog has a SUN at that position
        # [paralog, hg19_pos, ref, alt, hg38_pos]
        self.wl = {int(x[1]): [x[0], x[2], x[3], int(x[4])] for x in wl_list}

    def run(self):
        SunModel.run(self)


class IlpModel(object):
    def __init__(self, outDir, bpPenalty, dataPenalty, tightness, fastqFile, uuid, graph, localTempDir,
                 kmerSize, saveCounts=False, inferC=None, inferD=None):
        self.outDir = outDir
        self.uuid = uuid
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.fastqFile = fastqFile
        self.graph = graph
        self.localTempDir = localTempDir
        self.tightness = tightness
        self.saveCounts = saveCounts
        self.kmerSize = kmerSize
        self.inferC = inferC
        self.inferD = inferD

    def wigglePlots(self):
        explodedRawCounts = explodeResultDict(self.rawCounts, self.offsetMap)
        with open(os.path.join(self.outDir, "tracks", self.uuid + ".ILP.wig"), "w") as outf:
            outf.write(
                "track type=wiggle_0 name={} color=35,125,191 autoScale=off visibility=full alwaysZero=on yLineMark=2 "
                "viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n".format(self.uuid))
            for para in self.resultDict:
                for start, span, val in self.resultDict[para]:
                    outf.write("variableStep chrom=chr1 span={}\n".format(span))
                    outf.write("{} {}\n".format(start, val))
        with open(os.path.join(self.outDir, "tracks", self.uuid + ".KmerCounts.wig"), "w") as outf:
            outf.write(
                "track type=wiggle_0 name={} color=63,153,158 autoScale=off visibility=full alwaysZero=on yLineMark=2 "
                "viewLimits=0:4 yLineOnOff=on maxHeightPixels=100:75:50\n".format(self.uuid))
            for para in self.rawCounts:
                rawData = explodedRawCounts[para]
                windowedRawData = [sum(rawData[k:k+300])/300 for k in xrange(0, len(rawData), 300)]
                start = self.offsetMap[para]
                outf.write("fixedStep chrom=chr1 start={} step=300 span=300\n".format(start))
                for x in windowedRawData:
                    outf.write("{}\n".format(x))

    def run(self):
        if self.saveCounts is not True:
            countFile = os.path.join(self.localTempDir, self.uuid + ".Counts.fa")
        else:
            countFile = os.path.join(self.outDir, self.uuid + ".Counts.fa")
        if not os.path.exists(countFile):
            runJellyfish(self.localTempDir, countFile, self.fastqFile, self.uuid, self.kmerSize)
        G = pickle.load(open(self.graph, "rb"))
        dataCounts = Counter()
        normalizing = 0
        with open(countFile) as f:
            rm = ">\n"
            for count, seq in izip(*[f] * 2):
                seq = seq.translate(None, rm)
                count = int(count.translate(None, rm))
                rc = reverseComplement(seq)
                if seq in G.kmers or rc in G.kmers:
                    dataCounts[seq] += int(count)
                elif seq in G.normalizingKmers:
                    normalizing += int(count)
        normalizing /= (1.0 * len(G.normalizingKmers))
        P = KmerModel(G, normalizing, self.bpPenalty, self.dataPenalty, self.tightness, self.inferC, self.inferD)
        P.introduceData(dataCounts)
        P.solve()
        self.resultDict = P.reportCondensedCopyMap()
        self.rawCounts = P.reportCondensedNormalizedRawDataMap()
        self.offsetMap = P.getOffsets()
        self.wigglePlots()


def runJellyfish(localTempDir, countFile, fastqFile, uuid, kmerSize=49):
    """
    Runs jellyfish. -C flag is set to count both strands together.
    """
    jfFile = os.path.join(localTempDir, uuid + ".jf")
    system("jellyfish count -C -m {} -s 300M -o {} {}".format(kmerSize, jfFile, fastqFile))
    system("jellyfish dump -L 2 {} > {}".format(jfFile, countFile))


def downloadQuery(fastqPath, tempDir, key, queryString, uuid):
    """
    Downloads data from CGHub BAM Slicer
    """
    system(
        """curl --silent "{}" -u "{}" | samtools bamshuf -Ou /dev/stdin {} | samtools bam2fq /dev/stdin > {}""".format(
            queryString, "haussler:" + key, os.path.join(tempDir, "tmp"), fastqPath))
    if os.path.getsize(fastqPath) < 513:
        raise RuntimeError("curl did not download a BAM for {}. exiting.".format(uuid))


def alignQuery(fastqPath, remappedBamPath, tempDir, uuid, index):
    """
    Aligns to the notch locus
    """
    # align the extracted reads to the index
    sortedBamPath = os.path.join(tempDir, "{}.sorted".format(uuid))
    system("bwa mem -v 1 {} {} | samtools view -F 4 -bS - | samtools sort - {}".format(index, fastqPath, sortedBamPath))
    # samtools appends .bam to sorted bam files
    sortedBamPath += ".bam"
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


def combinedPlot(ilpDict, rawCounts, offsetMap, unfilteredSunDict, uuid, outDir):
    """
    Generates a final combined plot overlaying both ILP and SUN results.
    """
    colors = ["#9b59b6", "#3498db", "#e74c3c", "#34495e", "#2ecc71"]
    explodedRawCounts = explodeResultDict(rawCounts, offsetMap)
    explodedData = explodeResultDict(ilpDict, offsetMap)
    # used because the SUN model uses single letter labels
    paraMap = {"Notch2NL-A": "A", "Notch2NL-B": "B", "Notch2NL-C": "C", "Notch2NL-D": "D", "Notch2": "N"}
    sortedParalogs = ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D", "Notch2"]
    ax = plt.gca()
    fig, plots = plt.subplots(len(sortedParalogs), sharey=True)
    plt.yticks((0, 1, 2, 3, 4))
    plt.suptitle("kmer-DeBruijn ILP and SUN results")
    maxGap = max(len(x) for x in explodedRawCounts.itervalues())
    for i, (p, para) in enumerate(izip(plots, sortedParalogs)):
        data = explodedData[para]
        rawData = explodedRawCounts[para]
        windowedRawData = [sum(rawData[k:k+300])/300 for k in xrange(0, len(rawData), 300)]
        start = offsetMap[para]
        stop = start + len(data)
        rounded_max_gap = int(math.ceil(1.0 * maxGap / 10000) * 10000)
        p.axis([start, start + rounded_max_gap, 0, 4])
        x_ticks = range(start, start + rounded_max_gap, 20000)
        x_ticks.append(stop)
        p.axes.set_xticks(x_ticks)
        p.axes.set_xticklabels(map(str, range(0, (len(x_ticks) - 1) * 20000, 20000)) + [str(stop)])
        p.fill_between(range(start, stop), data, color=colors[i], alpha=0.8)
        p.plot(range(start, stop, 300), windowedRawData, alpha=0.8, linewidth=1.5)
        if len(unfilteredSunDict[paraMap[para]]) > 0:
            sunPos, sunVals = zip(*unfilteredSunDict[paraMap[para]])
            p.vlines(np.asarray(sunPos), np.zeros(len(sunPos)), sunVals, color="#E83535")
        for i in range(1, 4):
            p.axhline(y=i, ls="--", lw=0.8)
        p.set_title("{}".format(para))
    fig.subplots_adjust(hspace=0.8)
    plt.savefig(os.path.join(outDir, uuid + ".combined.png"), format="png")
    plt.close()

def explodeResultDict(rd, offsetMap):
    r = defaultdict(list)
    for para in rd:
        prevStart = offsetMap[para]
        for start, span, val in rd[para]:
            r[para].extend([val] * (start - prevStart))
            prevStart = start
    return r