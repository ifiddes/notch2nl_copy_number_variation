import sys, os, pysam, vcf, string, gzip
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
from lib.general_lib import formatRatio, rejectOutliers, opener

from jobTree.scriptTree.target import Target
from jobTree.src.bioio import system, logger, reverseComplement



class ModelWrapper(Target):
    """
    This Target runs all of the models.
    First, the fastq is extracted from the BAM slicer via curl and samtools.
    Next, the SUN model is ran (see SunModel)
    Then, the ILP model is ran (see IlpModel)
    Finally, the results of both models is used to build a combined plot.
    """
    def __init__(self, uuid, queryString, baseOutDir, bpPenalty, dataPenalty, keyFile, graph):
        Target.__init__(self)
        self.uuid = uuid
        self.queryString = queryString
        self.baseOutDir = baseOutDir
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.graph = graph
        self.key = open(keyFile).readline().rstrip()
        #index is a bwa index of the region to be aligned to (one copy of notch2)
        self.index = "./data/SUN_data/hs_n2.masked.fa"

    def downloadQuery(self):
        """
        Downloads data from CGHub BAM Slicer
        """
        fastqFile = os.path.join(self.getLocalTempDir(), self.uuid + ".fastq")
        logger.info("Downloading {} to {}".format(self.uuid, fastqFile))
        system("""curl --silent "{}" -u "{}" | samtools bamshuf -Ou /dev/stdin {} """
                """| samtools bam2fq /dev/stdin > {}""".format(self.queryString, 
                "haussler:" + self.key, os.path.join(self.getLocalTempDir(), "tmp"), fastqFile))
        if os.path.getsize(fastqFile) < 513:
            raise RuntimeError("curl did not download a BAM for {}. exiting.".format(self.uuid))
        #align the extracted reads to the index
        sortedBamPath = os.path.join(self.getLocalTempDir(), "{}.sorted".format(self.uuid))
        system("bwa mem -v 1 {} {} | samtools view -F 4 -bS - | samtools sort - {}".format(self.index, fastqFile, sortedBamPath))
        #samtools appends .bam to sorted bam files
        sortedBamPath += ".bam"
        #filter the SAM records and find the site coverage at each locus, creating VCFs
        remappedBamPath = os.path.join(self.getLocalTempDir(), "{}.remapped.sorted.bam".format(self.uuid))
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
        return remappedBamPath, fastqFile 

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
        plt.savefig(os.path.join(self.baseOutDir, self.uuid, self.uuid[:8] + ".combined.png"), format="png")
        plt.close()

    def run(self):
        bamPath, fastqFile = self.downloadQuery()
        sun = FilteredSunModel(self.baseOutDir, self.uuid, bamPath)
        sun.run()
        unfilteredSun = UnfilteredSunModel(self.baseOutDir, self.uuid, bamPath)
        unfilteredSun.run()
        ilp = IlpModel(self.baseOutDir, self.bpPenalty, self.dataPenalty, fastqFile, self.uuid, 
                self.graph, self.getLocalTempDir())
        ilp.run()
        self.combinedPlot(ilp.resultDict, sun.resultDict, unfilteredSun.resultDict, ilp.maxPos, ilp.offset)


class SunModel(object):
    """
    Runs the SUN model, where alelle fraction at unique sites is used to infer copy number
    """

    def findSiteCoverages(self, bamIn):
        """
        Runs mpileup on each SUN position and finds allele fractions
        """
        bases = set(["A", "T", "G", "C", "a", "t", "g", "c"])
        resultDict = defaultdict(list)
        nVals = []
        for pos, (para, ref, alt) in self.wl.iteritems():
            posStr = "chr1:{0}-{0}".format(pos)
            pileUp = pysam.mpileup("-q", "10", "-r", posStr, bamIn)
            if len(pileUp) == 0:
                continue
            pileUpStr = pileUp[0].split()
            if len(pileUpStr) != 6:
                continue
            pileUpResult = Counter(x.upper() for x in pileUpStr[4] if x in bases)
            if ref not in pileUpResult or alt not in pileUpResult:
                continue
            frac = formatRatio(pileUpResult[alt], sum(pileUpResult.values()))
            #invert fraction for Notch2 paralogs
            if para == "N":
                frac = 1 - frac
                nVals.append(frac)
            else:
                resultDict[para].append([pos, frac])
        norm = self.findNormalizingFactor(nVals)
        for para, result in resultDict.iteritems():
            resultDict[para] = [[x, 2.0 * y / norm] for x,y in result]
        return resultDict

    def findNormalizingFactor(self, r):
        """
        The values for Notch2NL-A-D are normalized by Notch2. Notch2 does not have CNV.
        """
        return np.mean(rejectOutliers(np.asarray(r)))

    def plotHistograms(self, resultDict):
        path = os.path.join(self.outDir, "{}.{}.png".format(self.uuid[:8], self.__class__.__name__))
        f, plts = plt.subplots(4, sharex=True)
        for para, p in izip(sorted(resultDict.keys()), plts):
            p.set_title("Notch2NL-{}".format(para))
            p.hist(np.asarray(zip(*resultDict[para])[1]), bins=40, range=(0.0,4.0), color='#1e90ff', normed=True)
            p.yaxis.tick_right()
            setp(p.get_yticklabels(), fontsize=9)
        plt.xticks((0, 1, 2, 3, 4))
        f.subplots_adjust(hspace=0.4)
        plt.savefig(path)

    def makeBedgraphs(self, resultDict):
        if not os.path.exists(os.path.join(self.outDir, "bedgraphs")):
            os.mkdir(os.path.join(self.outDir, "bedgraphs"))
        for para in resultDict:
            path = os.path.join(self.outDir, "bedgraphs", "{}.Notch2NL-{}.{}.bedGraph".format( 
                    self.uuid[:8], para, self.__class__.__name__))
            bedHeader = ("track type=bedGraph name={} autoScale=off visibility=full alwaysZero=on "
                    "yLineMark=0.2 viewLimits=0.0:0.4 yLineOnOff=on maxHeightPixels=100:75:50\n")
            with open(path, "w") as outf:
                outf.write(bedHeader.format(self.uuid[:8] + "_" + para))
                for pos, frac in resultDict[para]:
                    outf.write("\t".join(map(str, ["chr1", pos, pos + 1, frac])) + "\n")

    def run(self):
        self.resultDict = self.findSiteCoverages(self.bamPath)
        #plot the results
        self.plotHistograms(self.resultDict)
        self.makeBedgraphs(self.resultDict)
        #need to add (SUN-based) ILP here - hasn't been working with WGS data
        #self.call_ilp()
        #pickle.dump(self.resultDict, open(os.path.join(self.outDir, "resultDict.pickle"), "w"))


class UnfilteredSunModel(SunModel):
    def __init__(self, baseOutDir, uuid, bamPath):
        self.uuid = uuid
        self.bamPath = bamPath
        self.outDir = os.path.join(baseOutDir, self.uuid)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)
        #whitelist is a text file of whitelisted SUN positions - in this case, unfiltered
        with open("./data/SUN_data/unfiltered_whitelist.txt") as wl:
            wl_list = [x.split() for x in wl if not x.startswith("#")]
        #dict mapping genome positions to which paralog has a SUN at that position
        #[paralog, hg19_pos, ref, alt]
        self.wl = {int(x[1]) : [x[0], x[2], x[3]] for x in wl_list}

    def run(self):
        SunModel.run(self)


class FilteredSunModel(SunModel):
    def __init__(self, baseOutDir, uuid, bamPath):
        self.uuid = uuid
        self.bamPath = bamPath
        self.outDir = os.path.join(baseOutDir, self.uuid)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)
        #whitelist is a text file of whitelisted SUN positions - in this case, unfiltered
        with open("./data/SUN_data/whitelist.txt") as wl:
            wl_list = [x.split() for x in wl if not x.startswith("#")]
        #dict mapping genome positions to which paralog has a SUN at that position
        #[paralog, hg19_pos, ref, alt]
        self.wl = {int(x[1]) : [x[0], x[2], x[3]] for x in wl_list}

    def run(self):
        SunModel.run(self)


class IlpModel(object):
    def __init__(self, baseOutDir, bpPenalty, dataPenalty, fastqFile, uuid, graph, localTempDir):
        self.baseOutDir = baseOutDir
        self.uuid = uuid
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.fastqFile = fastqFile
        self.graph = graph
        self.localTempDir = localTempDir

        self.outDir = os.path.join(baseOutDir, self.uuid)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

    def plotResult(self, copyMap, maxPos, offset):
        #first do overlaid figure
        arrays = {}
        colors = ["#4D4D4D","#5DA5DA","#FAA43A","#60BD68","#F17CB0"]
        xvals = np.array(range(offset, offset + maxPos + 1),dtype="int")
        sortedParalogs = sorted(copyMap.keys())
        for para in sortedParalogs:
            arrays[para] = np.array(copyMap[para], dtype="int")
        fig = plt.figure()
        plt.axis([offset, offset + maxPos + 1, 0, 14])
        ax = plt.gca()
        total = sum(arrays.values())
        patches = []
        for i, para in enumerate(sortedParalogs):
            plt.fill_between(xvals, total, color=colors[i], alpha=0.7)
            total -= arrays[para]
            patches.append(mpatches.Patch(color=colors[i], label=para, alpha=0.7))
        plt.legend(patches, sortedParalogs)
        plt.suptitle("kmer-DeBruijn ILP results Notch2NL")
        plt.ylabel("Inferred Copy Number")
        plt.savefig(os.path.join(self.outDir, self.uuid[:8] + ".overlaid.ILP.png"), format="png")
        plt.close()
        #now do individual plots on one png
        fig, plots = plt.subplots(len(sortedParalogs), sharex=True, sharey=True)
        plt.suptitle("kmer-DeBruijn ILP results Notch2NL")
        plt.yticks((0, 1, 2, 3, 4))
        for i, (p, para) in enumerate(izip(plots, sortedParalogs)):
            p.axis([offset, offset + maxPos + 1, 0, 4])
            p.fill_between(xvals, copyMap[para], color=colors[i], alpha=0.7)
            p.set_title("{} Copy Number".format(para))
            p.scatter(xvals, copyMap[para], color=colors[i], alpha=0.7, s=0.3)
        fig.subplots_adjust(hspace=0.5)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        plt.savefig(os.path.join(self.outDir, self.uuid[:8] + ".separated.ILP.png"), format="png")
        plt.close()

    def run(self):
        jfFile = os.path.join(self.localTempDir, self.uuid + ".jf")
        countFile = os.path.join(self.localTempDir, self.uuid + ".Counts.fa")
        system("jellyfish count -m 49 -s 300M -o {} {}".format(jfFile, self.fastqFile))
        system("jellyfish dump -L 2 {} > {}".format(jfFile, countFile))
        
        G = pickle.load(open(self.graph, "rb"))
        with opener(countFile) as f:
            #using string translate has been shown to be faster than using any other means
            #of removing characters from a string
            rm = ">\n"
            dataCounts = {}
            for count, seq in izip(*[f]*2):
                seq = seq.translate(None, rm)
                rc = reverseComplement(seq)
                if seq in G.kmers or seq in G.normalizingKmers:
                    dataCounts[seq] = int(count.translate(None, rm))
                elif rc in G.reverseKmers or rc in G.reverseNormalizingKmers:
                    dataCounts[rc] = int(count.translate(None, rm))

        #adjust ILP penalties for coverage in this sequencing run
        normalizing = (( 1.0 * sum(dataCounts.get(x, 0) for x in G.normalizingKmers) 
                + sum(dataCounts.get(x, 0) for x in G.reverseNormalizingKmers) ) 
                / len(G.normalizingKmers))

        P = KmerModel(G, normalizing, self.bpPenalty, self.dataPenalty)
        P.introduceData(dataCounts)
        P.solve()
        self.resultDict, self.maxPos = P.reportCopyMap()
        self.plotResult(self.resultDict, self.maxPos, G.offset)
        self.offset = G.offset
