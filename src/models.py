import sys, os, pysam, vcf, string
import cPickle as pickle
from itertools import izip
from collections import defaultdict, Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import setp

from src.kmerModel import KmerModel
from lib.general_lib import formatRatio, rejectOutliers

from jobTree.scriptTree.target import Target
from jobTree.src.bioio import system, logger, reverseComplement


class ModelWrapper(Target):
    def __init__(self, uuid, queryString, baseOutDir, bpPenalty, dataPenalty, keyFile, graph):
        Target.__init__(self)
        self.uuid = uuid
        self.queryString = queryString
        self.baseOutDir = baseOutDir
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.graph = graph
        self.key = open(keyFile).readline().rstrip()

    def downloadQuery(self):
        fastqFile = os.path.join(self.getLocalTempDir(), self.uuid + ".fastq")
        logger.info("Downloading {} to {}".format(self.uuid, fastqFile))
        system("""curl --silent "{}" -u "{}" | samtools bamshuf -Ou /dev/stdin {} """
                """| samtools bam2fq /dev/stdin > {}""".format(self.queryString, 
                "haussler:" + self.key, os.path.join(self.getLocalTempDir(), "tmp"), fastqFile))
        return fastqFile

    def combinedPlot(self, ilpDict, sunDict, maxPos, offset, normalizing):
        arrays = {}
        colors = ["#4D4D4D","#5DA5DA","#FAA43A","#60BD68","#F17CB0"]
        xvals = np.array(range(offset, offset + maxPos + 1),dtype="int")
        sortedParalogs = sorted(ilpDict.keys())
        for para in sortedParalogs:
            arrays[para] = np.array(ilpDict[para], dtype="int")
        ax = plt.gca()
        fig, plots = plt.subplots(len(sortedParalogs), sharex=True, sharey=True)
        plt.yticks((0, 1, 2, 3))
        for i, (p, para) in enumerate(izip(plots, sortedParalogs)):
            p.axis([offset, offset + maxPos + 1, 0, 3])
            p.fill_between(xvals, ilpDict[para], color=colors[i], alpha=0.7)
            p.set_title("{} Copy Number".format(para))
            sunPos, sunVals = zip(*sunDict[para])
            plt.vlines(np.asarray(sunPos), np.zeros(len(sunPos)), np.asarray(sunVals), lw=2)
        fig.subplots_adjust(hspace=0.5)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False) 
        plt.savefig(os.path.join(self.baseOutDir, self.uuid + ".sun.ilp.png"), format="png")
        plt.close()     

    def run(self):
        fastqFile = self.downloadQuery()
        if os.path.getsize(fastqFile) < 513:
            raise RuntimeError("curl did not download a BAM for {}. exiting.".format(self.uuid))
        sun = FilteredSunModel(self.baseOutDir, fastqFile, self.uuid, self.getLocalTempDir())
        sunDict, normalizing = sun.run()
        ilp = IlpModel(self.baseOutDir, self.bpPenalty, self.dataPenalty, fastqFile, self.uuid, 
                self.graph, self.getLocalTempDir())
        ilpDict, maxPos, offset = ilp.run()
        self.combinedPlot(ilpDict, sunDict, maxPos, offset)


class SunModel(object):
    def filterAndIndex(self, bamIn, bamOut):
        header = {"HD": {"VN": "1.3"}, "SQ": [{"LN": 248956422, "SN": "chr1"}]}

        outfile = pysam.Samfile(bamOut, "wb", header=header)
        bamfile = pysam.Samfile(bamIn, "rb")  

        for record in bamfile:
            if not record.is_unmapped:
                chrom, span = bamfile.getrname(record.tid).split(":")
                start, end = map(int, span.split("-"))
                record.pos = record.pos + start - 1
                outfile.write(record)

        outfile.close()
        system("samtools index {}".format(bamOut))

    def findSiteCoverages(self, bamIn):
        bases = set(["A", "T", "G", "C", "a", "t", "g", "c"])
        resultDict = defaultdict(list)
        
        for pos, (para, ref, alt) in self.wl.iteritems():
            posStr = "chr1:{0}-{0}".format(pos)
            pileUpStr = pysam.Mpileup("-q", "10", "-r", posStr, bamIn)[0].split()
            if len(pileUpStr) != 6:
                continue

            pileUpResult = Counter(x.upper() for x in pileUpStr[4] if x in bases)
            if ref not in pileUpResult or alt not in pileUpResult:
                continue

            frac = formatRatio(pileUpResult[alt], sum(pileUpResult.values()))
            #invert fraction for Notch2 paralogs
            if para == "N":
                frac = 1 - frac
            resultDict[para].append([pos, frac])
        return resultDict

    def findNormalizingFactor(self, r):
        return np.mean(rejectOutliers(np.asarray(r))) / 0.8

    def plotHistograms(self, resultDict, path):
        f, plts = plt.subplots(5, sharex=True)
        fig = plt.figure()
        for para, p in izip(sorted(resultDict.keys()), plts):
            v = np.asarray(zip(*resultDict[para])[1]) * self.normalizing
            p.set_ylabel("{}".format(para), size=9)
            p.hist(v, bins=30, range=(0.0,1.0), color='#1e90ff', normed=True)
            p.yaxis.tick_right()
            setp(p.get_yticklabels(), fontsize=8)
        plt.savefig(path)

    def makeBedgraphs(self, resultDict):
        for para in resultDict:
            path = os.path.join(self.outDir,"bedgraphs", "{}.{}.bedGraph".format(self.uuid, 
                    self.header))
            bedHeader = ("track type=bedGraph name={} autoScale=off visibility=full alwaysZero=on "
                    "yLineMark=0.2 viewLimits=0.0:0.4 yLineOnOff=on maxHeightPixels=100:75:50\n")
            with open(path, "w") as outf:
                outf.write(bedHeader.format(name[:8] + "_" + para))
                for pos, frac in resultDict[para]:
                    outf.write("\t".join(map(str, ["chr1", pos, pos + 1, frac * self.normalizing])) \
                            + "\n")

    def run(self):
        #align the extracted reads to the index
        sortedBamPath = os.path.join(self.localTempDir, "{}.{}.sorted".format(self.uuid, self.header))
        system("bwa mem -v 1 {} {} | samtools view -bS - | samtools sort - {}".format(self.index, 
                self.fastqFile, sortedBamPath))
        #samtools appends .bam to sorted bam files
        sortedBamPath += ".bam"
        
        #filter the SAM records and find the site coverage at each locus, creating VCFs
        remappedBamPath = os.path.join(self.localTempDir, 
                "{}.{}.remapped.sorted.bam".format(self.uuid, self.header))
        resultDict = self.findSiteCoverages(remappedBamPath)     
        self.normalizing = self.findNormalizingFactor(zip(*resultDict["N"])[1])
        #plot the results
        self.plotHistograms(resultDict, os.path.join(self.outDir, "{}.{}.png".format(self.uuid, 
                self.header)))
        self.makeBedgraphs(resultDict)
        #need to add (SUN-based) ILP here - hasn't been working with WGS data
        #self.call_ilp()
        return resultDict, self.normalizing


class UnfilteredSuns(SunModel):
    def __init__(self, baseOutDir, fastqFile, uuid, localTempDir):
        self.uuid = uuid
        self.fastqFile = fastqFile
        self.localTempDir = localTempDir

        self.outDir = os.path.join(baseOutDir, self.uuid)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

        ###########
        # these are hard coded files. Shitty but whatever.
        ###########

        self.header = "unfilteredSun"
        #index is a bwa index of the region to be aligned to (one copy of notch2)
        self.index = "./data/index/hs_n2.masked.fa"

        #whitelist is a text file of whitelisted SUN positions - in this case, unfiltered
        with open("./data/unfiltered_whitelist.txt") as wl:
            wl_list = [x.split() for x in wl if not x.startswith("#")]
        #dict mapping genome positions to which paralog has a SUN at that position
        #[paralog, hg19_pos, ref, alt]
        self.wl = {int(x[1]) : [x[0], x[2], x[3]] for x in wl_list}

    def run(self):
        SunModel.run(self)


class FilteredSuns(SunModel):
    def __init__(self, baseOutDir, fastqFile, uuid, localTempDir):
        self.uuid = uuid
        self.fastqFile = fastqFile
        self.localTempDir = localTempDir

        self.outDir = os.path.join(baseOutDir, self.uuid)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

        ###########
        # these are hard coded files. Shitty but whatever.
        ###########

        self.header = "filteredSun"
        #index is a bwa index of the region to be aligned to (one copy of notch2)
        self.index = "./data/index/hs_n2.masked.fa"

        #whitelist is a text file of whitelisted SUN positions - in this case, unfiltered
        with open("./data/whitelist.txt") as wl:
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
        plt.savefig(os.path.join(self.baseOutDir, self.uuid + ".png"), format="png")
        plt.close()
        #now do individual plots on one png
        fig, plots = plt.subplots(len(sortedParalogs), sharex=True, sharey=True)
        plt.yticks((0, 1, 2, 3))
        for i, (p, para) in enumerate(izip(plots, sortedParalogs)):
            p.axis([offset, offset + maxPos + 1, 0, 3])
            p.fill_between(xvals, copyMap[para], color=colors[i], alpha=0.7)
            p.set_title("{} Copy Number".format(para))
            p.scatter(xvals, copyMap[para], color=colors[i], alpha=0.7, s=0.3)
        fig.subplots_adjust(hspace=0.5)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        plt.savefig(os.path.join(self.baseOutDir, self.uuid + ".separated.png"), format="png")
        plt.close()

    def run(self):
        jfFile = os.path.join(self.localTempDir, self.uuid + ".jf")
        countFile = os.path.join(self.localTempDir, self.uuid + ".Counts.fa")
        system("jellyfish count -m 49 -s 300M -o {} {}".format(jfFile, self.fastqFile))
        system("jellyfish dump -L 2 {} > {}".format(jfFile, countFile))
        
        G = pickle.load(self.graph)
        with open(countFile) as f:
            #using string translate has been shown to be faster than using any other means
            #of removing characters from a string
            rm = ">\n"
            dataCounts = {}
            for count, seq in izip(*[f]*2):
                seq = seq.translate(None, rm)
                rc = reverseComplement(seq)
                if seq in G.kmers or seq in G.normalizing_kmers:
                    dataCounts[seq] = int(count.translate(None, rm))
                elif rc in G.reverse_kmers or rc in G.reverse_normalizing_kmers:
                    dataCounts[rc] = int(count.translate(None, rm))

        #adjust ILP penalties for coverage in this sequencing run
        normalizing = (( 1.0 * sum(dataCounts.get(x, 0) for x in G.normalizing_kmers) 
                + sum(dataCounts.get(x, 0) for x in G.reverse_normalizing_kmers) ) 
                / len(G.normalizing_kmers))

        P = KmerModel(G, normalizing, self.bpPenalty, self.dataPenalty)
        P.introduceData(dataCounts)
        P.solve()
        copyMap, maxPos = P.reportCopyMap()
        plotResult(copyMap, maxPos, G.offset)
        return copyMap, maxPos, G.offset
