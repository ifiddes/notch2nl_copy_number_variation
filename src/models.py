import sys, os, pysam, vcf, string
import cPickle as pickle
from itertools import izip

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import setp

from src.kmerModel import KmerModel
from src.sunIlpModel import sunIlpModel

from jobTree.scriptTree.target import Target
from jobTree.src.bioio import system, reverseComplement


class ModelWrapper(Target):
    def __init__(self, uuid, queryString, baseOutDir, bpPenalty, dataPenalty, fullSun, originalSun, 
                Ilp, keyFile):
        Target.__init__(self)
        self.uuid = uuid
        self.queryString = queryString
        self.baseOutDir = baseOutDir
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.doFullSun = fullSun
        self.doOriginalSun = originalSun
        self.doIlp = Ilp
        self.key = open(keyFile).readline().rstrip()

    def downloadQuery(self):
        self.fastqFile = os.path.join(self.getLocalTempDir(), self.uuid + ".fastq")
        system(("curl --silent {queryString} -u {key} | samtools bamshuf -Ou /dev/stdin {tmpDir}/tmp"
                " | samtools bam2fq /dev/stdin > {fastqFile}").format(key=self.key, 
                queryString=self.queryString, fastqFile=self.fastqFile))

    def run(self):
        self.downloadQuery()
        if self.doFullSun is True:
            Target.addChildTarget(FullSunModel(self.baseOutDir, self.fastqFile, self.uuid))
        if self.doOriginalSun is True:
            Target.addChildTarget(OriginalSunModel(self.baseOutDir, self.fastqFile, self.uuid))
        if self.doIlp is True:
            Target.addChildTarget(IlpModel(self.baseOutDir, self.bpPenalty, self.dataPenalty, 
                    self.fastqFile, self.uuid))


class FullSunModel(AbstractSunModel):
    def __init__(self, baseOutDir, fastqFile, uuid):
        Target.__init__(self)

        self.uuid = uuid
        self.fastqFile = fastqFile

        self.outDir = os.path.join(baseOutDir, self.uuid)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

        ###########
        # these are hard coded files. Shitty but whatever.
        ###########

        #his is used to differentiate output files from the OriginalSunModel
        self.header = "fullSun"
        #index is a bwa index of the region to be aligned to (one copy)
        self.index = "./data/index/hg38_n2.masked.fa"
        #whitelist is a text file of whitelisted SUN positions for the FullSunModel
        with open("./data/textfiles/full_whitelist.txt") as wl:
            wl_list = [x.split() for x in wl if not x.startswith("#")]
        #dict mapping genome positions to which paralog has a SUN at that position
        self.wl = {x[1] : x[0] for x in wl}
        #vcf files used for this model
        self.vcfs = {x : os.path.join("data", "full_vcfs", x + ".vcf.gz") \
                for x in ["A", "B", "C", "D", "N"]}
        #full chr1 fasta sequence used for this model
        self.chr1 = "./data/fastas/hg38_chr1_repeat_masked.fa"

    def run(self):
        AbstractSunModel.run(self)

class OriginalSunModel(AbstractSunModel):
    def __init__(self, baseOutDir, fastqFile, uuid):
        Target.__init__(self)
        self.uuid = uuid
        self.fastqFile = fastqFile

        self.outDir = os.path.join(baseOutDir, self.uuid)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

        ###########
        # these are hard coded files. Shitty but whatever.
        ###########

        #his is used to differentiate output files from the OriginalSunModel
        self.header = "originalSun"
        #index is a bwa index of the region to be aligned to (one copy)
        self.index = "./data/index/hs_n2.masked.fa"
        #whitelist is a text file of whitelisted SUN positions for the FullSunModel
        with open("./data/textfiles/original_whitelist.txt") as wl:
            wl_list = [x.split() for x in wl if not x.startswith("#")]
        #dict mapping genome positions to which paralog has a SUN at that position
        self.wl = {x[1] : x[0] for x in wl_list}
        #vcf files used for this model
        self.vcfs = {x : os.path.join("data", "original_vcfs", x + ".vcf.gz") \
                for x in ["A", "B", "C", "D", "N"]}
        #full chr1 fasta sequence used for this model
        self.chr1 = "./data/fastas/hg19_ucsc.fa"

    def run(self):
        AbstractSunModel.run(self)


class AbstractSunModel(Target):

    def filterSam(self, bamIn, bamOut):
        header = {"HD": {"VN": "1.3"}, "SQ": [{"LN": 248956422, "SN": "chr1"}]}

        outfile = pysam.Samfile(bamOut, "wb", header=header)
        samfile = pysam.Samfile(bamIn, "rb")  

        for record in samfile:
            if not record.is_unmapped:
                chrom, span = samfile.getrname(record.tid).split(":")
                start, end = map(int, span.split("-"))
                record.pos = record.pos + start - 1
                outfile.write(record)

        outfile.close()

    def runGetSiteCoverage(self, bamIn):
        #I am not refactoring this legacy code, sorry.
        outVcfs = []
        for para, vcf in self.vcfs.iteritems():
            outVcfs.append(os.path.join(self.getLocalTempDir(), "{uuid}.{header}.{para}.vcf".format(
                    uuid=self.uuid, header=self.header, para=para)))
            system(("python ./src/getsitecoverage.py -t SNV -b {in} -v {vcf} "
                    "-i {chr1} --chr > {out}").format(in=bamIn, vcf=vcf, chr1=self.chr1,
                    out=outVcfs[-1]))
        return outVcfs

    def plot_histograms(self,data, path):
        f, plts = plt.subplots(5, sharex=True)
        fig = plt.figure()
        for para, p in izip(sorted(data.keys()), plts):
            v = zip(*data[para])[1]
            p.set_ylabel("{}".format(para), size=9)
            p.hist(1-np.asarray(v), bins=30, range=(0.0,1.0), color='#1e90ff', normed=True)
            p.yaxis.tick_right()
            setp(p.get_yticklabels(), fontsize=8)
        plt.savefig(path)

    def makeHg38Bedgraphs(self, resultDict):
        for para in resultDict:
            path = os.path.join(self.outDir, "{uuid}.{header}.bedGraph".format(uuid=self.uuid,
                    header=self.header))
            bedHeader = ("track type=bedGraph name={} autoScale=off visibility=full alwaysZero=on "
                    "yLineMark=0.2 viewLimits=0.0:0.4 yLineOnOff=on maxHeightPixels=100:75:50\n")
            with open(path, "w") as outf:
                outf.write(bedHeader.format(name[:6] + "_" + para))
                for hg38_pos, frac in resultDict[para]:
                    outf.write("\t".join(map(str, ["chr1", hg38_pos-1, hg38_pos, frac]))+"\n")

    def run(self):
        #align the extracted reads to the index
        sortedBamPath = os.path.join(self.getLocalTempDir(), "{uuid}.{header}.sorted.bam".format(
                    uuid=self.uuid, header=self.header))
        system("bwa mem -v 1 {index} {fastq} | samtools view -bS - | samtools sort - {out}".format(
                    index=self.index, fastq=self.fastqFile, out=sortedBamPath))
        #filter the SAM records and find the site coverage at each locus, creating VCFs
        remappedBamPath = os.path.join(self.getLocalTempDir(), 
                "{uuid}.{header}.remapped.sorted.bam".format(uuid=self.uuid, header=self.header))
        self.filterSam(sortedBamPath, remappedBamPath)
        outVcfs = self.runGetSiteCoverage(remappedBamPath)        

        #loop over these created VCFs and find allele fraction of SUNs
        resultDict = {"A":[], "B":[], "C":[], "D":[], "N":[]}
        for v in outVcfs:
            v = vcf.Reader(file(v))
            for record in v:
                if str(record.POS) in self.wl:
                    paralog = wl[str(record.POS)]
                    if paralog in resultDict:
                        val = float(record.INFO["ALTFRAC"][0])
                        resultDict[paralog].append((record.POS, val))

        self.plot_histograms(resultDict, os.path.join(self.outDir, 
                "{uuid}.{header}.png".format(uuid=self.uuid, header=self.header)))
        self.makeHg38Bedgraphs(resultDict)
        #need to add (SUN-based) ILP here - hasn't been working with WGS data
        #self.call_ilp()


class IlpModel(Target):
    def __init__(self, baseOutDir, bpPenalty, dataPenalty, fastqFile, uuid, graph):
        Target.__init__(self)
        self.baseOutDir = baseOutDir
        self.uuid = uuid
        self.bpPenalty = bpPenalty
        self.dataPenalty = dataPenalty
        self.fastqFile = fastqFile
        self.graph = graph

        self.outDir = os.path.join(baseOutDir, self.uuid)
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)

    def plot_result(self, copy_map, max_pos):
        #first do overlaid figure
        arrays = {}
        colors = ["#4D4D4D","#5DA5DA","#FAA43A","#60BD68","#F17CB0"]
        xvals = np.array(range(max_pos+1),dtype="int")
        sorted_paralogs = sorted(copy_map.keys())
        for para in sorted_paralogs:
            arrays[para] = np.array(copy_map[para], dtype="int")
        fig = plt.figure()
        plt.axis([0, max_pos+1, 0, 14])
        ax = plt.gca()
        total = sum(arrays.values())
        patches = []
        for i, para in enumerate(sorted_paralogs):
            plt.fill_between(xvals, total, color=colors[i], alpha=0.7)
            total -= arrays[para]
            patches.append(mpatches.Patch(color=colors[i], label=para, alpha=0.7))
        plt.legend(patches, sorted_paralogs)
        plt.suptitle("kmer-DeBruijn ILP results Notch2NL")
        plt.ylabel("Inferred Copy Number")
        plt.savefig(os.path.join(self.baseOutDir, self.uuid + ".png"), format="png")
        plt.close()
        #now do individual plots on one png
        fig, plots = plt.subplots(len(sorted_paralogs), sharex=True, sharey=True)
        plt.yticks((0, 1, 2, 3))
        for i, (p, para) in enumerate(izip(plots, sorted_paralogs)):
            p.axis([0, max_pos+1, 0, 3])
            p.fill_between(xvals, copy_map[para], color=colors[i], alpha=0.7)
            p.set_title("{} Copy Number".format(para))
            p.scatter(xvals, copy_map[para], color=colors[i], alpha=0.7, s=0.3)
        fig.subplots_adjust(hspace=0.5)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        plt.savefig(os.path.join(self.baseOutDir, self.uuid + ".separated.png"), format="png")
        plt.close()

    def run(self):
        jfFile = os.path.join(self.getGLobalTempDir(), self.uuid + ".jf")
        countFile = os.path.join(self.getLocalTempDir(), self.uuid + ".Counts.fa")
        system("jellyfish count -m 49 -s 300M -o {jfFile} {fastq}".format(jfFile=jfFile,
                fastqFile=self.fastqFile))
        system("jellyfish dump -L 2 {jfFile} > {countFile}".format(jfFile=jfFile, 
                countFile=countFile))
        
        G = pickle.load(self.graph)
        with open(countFile) as f:
            #using string translate has been shown to be faster than using any other means
            #of removing characters from a string
            rm = ">\n"
            data_counts = {}
            for count, seq in izip(*[f]*2):
                seq = seq.translate(None, rm)
                rc = reverseComplement(seq)
                if seq in G.kmers or seq in G.normalizing_kmers:
                    data_counts[seq] = int(count.translate(None, rm))
                elif rc in G.reverse_kmers or rc in G.reverse_normalizing_kmers:
                    data_counts[rc] = int(count.translate(None, rm))

        #adjust ILP penalties for coverage in this sequencing run
        normalizing = 1.0 * (( sum(data_counts.get(x, 0) for x in G.normalizing_kmers) 
                + sum(data_counts.get(x, 0) for x in G.reverse_normalizing_kmers) ) 
                / len(G.normalizing_kmers))

        P = KmerModel(G.paralogs, normalizing, self.bpPenalty, self.dataPenalty)
        P.build_blocks(G)
        P.introduce_data(data_counts)
        P.solve()
        copy_map, max_pos = P.report_copy_map()
        plot_result(copy_map, max_pos)
