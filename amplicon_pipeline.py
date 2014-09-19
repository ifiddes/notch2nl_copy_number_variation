#!/usr/bin/env python

"""Re-written python-only version of the notch2nl CNV inference pipeline.
Designed for use with amplicon sequencing."""

import sys, os, subprocess, argparse, pysam, vcf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import setp
from itertools import izip
from collections import OrderedDict as od

header = "track type=bedGraph name={} autoScale=off visibility=full alwaysZero=on yLineMark=0.2 viewLimits=0.0:0.4 yLineOnOff=on maxHeightPixels=100:75:50\n"

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--fwd", type=str, help="forward read", required=True)
    parser.add_argument("--rev", type=str, help="reverse read", required=True)
    parser.add_argument("--out", type=str, help="output folder", default="remapped")
    parser.add_argument("--name", type=str, help="experiment name", required=True)
    parser.add_argument("--whitelist", type=argparse.FileType("r"), help="whitelist file", default="whitelist.txt")
    parser.add_argument("--cores", type=int, help="# of cores for alignments", default=10)
    parser.add_argument("--index", type=str, help="base (ends in .fa) of bwa index of n2 locus", default="index/hs_n2.masked.fa")
    return parser.parse_args()


def call_bwa_samtools(fwd, rev, basename, cores, index):
    call = map(str, ["sh", "run_bwa_samtools.sh", fwd, rev, basename, cores, index])
    subprocess.Popen(call).wait()

def call_ilp(basename, name):
    Avcf, Bvcf, Cvcf, ABvcf = [os.path.join(basename, name + x + ".vcf") for x in ['.A','.B','.C','.AB']]
    pngpath = os.path.join(basename, name + "_ILP.png")
    call = ["python", "integer_linear_programming_notch2nl.py", "--A", Avcf, "--B", Bvcf, "--C", Cvcf,
            "--AB", ABvcf, "--png", pngpath, "--name", name, "--normalize", "--save_lp_results"]
    subprocess.Popen(call).wait()

def relabel(ax):
    labels = []
    for tick in ax.get_yticks():
        labels.append(float(tick))
    labels = map(str, np.asarray(labels)/100.0)
    return labels


def plot_histograms(data, path):
    N = zip(*data['N'])[1]
    A = zip(*data['A'])[1]
    B = zip(*data['B'])[1]
    C = zip(*data['C'])[1]
    D = zip(*data['D'])[1]
    fig = plt.figure()
    f, (axN2, axA, axB, axC, axD) = plt.subplots(5, sharex=True)
    axN2.set_ylabel("NOTCH2", size=9)
    axN2.hist(1-np.asarray(N), bins=30, range=(0.0,1.0), color='#1e90ff', normed=True)
    axN2.yaxis.tick_right()
    setp(axN2.get_yticklabels(), fontsize=8)
    axN2.set_yticklabels(relabel(axN2))
    axA.set_ylabel("NOTCH2NL-A", size=9)
    axA.hist(np.asarray(A), bins=30, range=(0.0,1.0), color='#1e90ff', normed=True)
    axA.yaxis.tick_right()
    setp(axA.get_yticklabels(), fontsize=8)
    axA.set_yticklabels(relabel(axA))
    axB.set_ylabel("NOTCH2NL-B", size=9)
    axB.hist(np.asarray(B), bins=30, range=(0.0,1.0), color='#1e90ff', normed=True)
    axB.yaxis.tick_right()
    setp(axB.get_yticklabels(), fontsize=8)
    axB.set_yticklabels(relabel(axB))
    axC.set_ylabel("NOTCH2NL-C", size=9)
    axC.hist(np.asarray(C), bins=30, range=(0.0,1.0), color='#1e90ff', normed=True)
    axC.yaxis.tick_right()
    setp(axC.get_yticklabels(), fontsize=8)
    axC.set_yticklabels(relabel(axC))
    axD.set_ylabel("NOTCH2NL-D", size=9)
    axD.hist(np.asarray(D), bins=30, range=(0.0,1.0), color='#1e90ff', normed=True)
    axD.yaxis.tick_right()
    setp(axD.get_yticklabels(), fontsize=8)
    axD.set_yticklabels(relabel(axD))
    plt.savefig(path)

def make_bedgraph(result_dict, wl, path, name):
    outf = open(path, "w")
    outf.write(header.format(name))
    for para in result_dict:
        for hg19_pos, frac in result_dict[para]:
            if para == "N":
                frac = 1 - frac
            elif para == "AB":
                frac = frac / 2
            chm1_pos = wl[str(hg19_pos)][2]
            outf.write("\t".join(map(str, ["chr1", int(chm1_pos)-1, chm1_pos, frac]))+"\n")
    outf.close()



def main(args):
    args = parse_args(args)
    #base name used for lots of things
    if not os.path.isdir(os.path.join(args.out, args.name)):
        os.mkdir(os.path.join(args.out, args.name))

    basename = os.path.join(args.out, args.name)

    #call bwa and run samtools mpileup
    call_bwa_samtools(args.fwd, args.rev, os.path.join(basename, args.name), args.cores, args.index)
    
    #make set of wl positions
    wl = [x.split() for x in args.whitelist][1:]
    wl = {x[0]:(x[1],x[2],x[3]) for x in wl}

    vcf_paths = [os.path.join(basename, args.name + x + ".vcf") for x in ['.N','.A','.B','.C','.D','.AB']]

    result_dict = {"A":list(), "B":list(), "C":list(), "D":list(), "N":list(), "AB":list()}
    for v in vcf_paths:
        v = vcf.Reader(file(v))
        for record in v:
            if str(record.POS) in wl:
                paralog, weight, chm1_pos = wl[str(record.POS)]
                result_dict[paralog].append((record.POS, float(weight) * float(record.INFO["ALTFRAC"][0])))

    plot_histograms(result_dict, os.path.join(basename, args.name + "_hist.png"))
    make_bedgraph(result_dict, wl, os.path.join(basename, args.name + ".bedGraph"), args.name)
    call_ilp(basename, args.name)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
