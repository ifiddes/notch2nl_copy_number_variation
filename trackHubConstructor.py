#!/usr/bin/env python2.7

"""

Digs through the output folder specified and builds a trackHub out of it in the target directory.

"""

import sys, os, argparse
from jobTree.src.bioio import system
from lib.general_lib import DirType

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "-o", type=DirType, default="./output/",
            help="base output directory that will be read. Default is ./output/")
    parser.add_argument("--assembly_dir", "-a", type=DirType, required=True,
            help="Dir to put assembly hub into.")
    parser.add_argument("--name", required=True,
            help="assembly hub name")
    parser.add_argument("--hg38_chrom_sizes", help="hg38 chrom.sizes file", default="hg38.chrom.sizes")
    parser.add_argument("--hg19_chrom_sizes", help="hg19 chrom.sizes file", default="hg19.chrom.sizes")
    return parser.parse_args()


def startHub(d, name):
    with open(os.path.join(d, "hub.txt"), "w") as outf:
        outf.write("hub {}\nshortLabel {}\nlongLabel {}\ngenomesFile genomes.txt\nemail ian.t.fiddes@gmail.com\n".format(name, name + " Notch2NL", name + " Notch2NL"))
    with open(os.path.join(d, "genomes.txt"), "w") as outf:
        outf.write("genome hg38\ntrackDb hg38/trackDb.txt\ndefaultPos chr1:145987622-149723055\n")
        outf.write("genome hg19\ntrackDb hg19/trackDb.txt\ndefaultPos chr1:120414641-120651853\n")
    if not os.path.exists(os.path.join(d, "hg38")):
        os.mkdir(os.path.join(d, "hg38"))
    if not os.path.exists(os.path.join(d, "hg19")):
        os.mkdir(os.path.join(d, "hg19"))
    with open(os.path.join(d, "hg38", "Notch2NL.html"), "w") as outf:
         outf.write("Notch2NL {}\n".format(name))
    with open(os.path.join(d, "hg19", "Notch2NL.html"), "w") as outf:
         outf.write("Notch2NL {}\n".format(name))        

def buildTrackDb(d, paths):
    with open(os.path.join(d, "hg38", "trackDb.txt"), "w") as outf:
        for g, [sun, ilp, sun_ilp, ilp_counts, hg19_sun_paths] in paths.iteritems():
            if ilp is not None:
            #    outf.write("track {0}_combined\ncontainer multiWig\nshortLabel {0} Combined\nlongLabel {0} Combined\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\n\n".format(g))
            #    outf.write("track {0}_SUN_combined\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_combined\n\n".format(g, os.path.basename(sun)))
            #    outf.write("track {0}_ILP_combined\ncolor 242,148,176\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_combined\n\n\n".format(g, os.path.basename(ilp)))
            #    outf.write("track {0}_raw_counts_combined\ncolor 122,113,116\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_kmer_ILP\n\n".format(g, os.path.basename(ilp_counts)))                
                if sun_ilp is not None:
                    outf.write("track {0}_SUN_ILP_combined\ncolor 191,191,191\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_combined\n\n\n".format(g, os.path.basename(sun_ilp)))
                outf.write("track {0}_kmer_ILP\ncontainer multiWig\nshortLabel {0} k-mer ILP\nlongLabel {0} k-mer ILP\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\n\n".format(g))
                outf.write("track {0}_SUN_kmer\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_kmer_ILP\n\n".format(g, os.path.basename(sun)))
                outf.write("track {0}_ILP_kmer\ncolor 242,148,176\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_kmer_ILP\n\n\n".format(g, os.path.basename(ilp)))
                outf.write("track {0}_raw_counts\ncolor 137,137,137\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_kmer_ILP\n\n".format(g, os.path.basename(ilp_counts)))
            if sun_ilp is not None and sun is not None:
                outf.write("track {0}_SUN_ILP\ncontainer multiWig\nshortLabel {0} SUN ILP\nlongLabel {0} SUN ILP\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\n\n".format(g))
                outf.write("track {0}_SUN_sun\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_SUN_ILP\n\n".format(g, os.path.basename(sun)))
                outf.write("track {0}_SUN_ILP_sun\ncolor 191,191,191\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_SUN_ILP\n\n\n".format(g, os.path.basename(sun_ilp)))
    with open(os.path.join(d, "hg19", "trackDb.txt"), "w") as outf:
        for g, [sun, ilp, sun_ilp, ilp_counts, hg19_sun_paths] in paths.iteritems():
            outf.write("track {0}_SUN\ncontainer multiWig\nshortLabel {0} SUN\nlongLabel {0} SUN\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\n\n".format(g))
            for p in hg19_sun_paths:
                outf.write("track {0}_SUN_{2}\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}_SUN_ILP\n\n".format(g, os.path.basename(p), os.path.basename(p).split(".")[0]))


def main():
    args = parse_args()
    startHub(args.assembly_dir, args.name)

    genomes = [x for x in os.listdir(args.output)]
    paths = {}
    for g in genomes:
        bg = os.path.join(args.output, g, "tracks", g + ".UnfilteredSunModel.hg38.bedGraph")
        ilp = os.path.join(args.output, g, "tracks", g + ".ILP.wig")
        sun = os.path.join(args.output, g, "tracks", g + ".UnfilteredSunModel.hg38.SUN_ILP.wiggle")
        ilp_counts = os.path.join(args.output, g, "tracks", g + ".KmerCounts.wig")

        out_sun_bw = os.path.join(args.assembly_dir, "hg38", g + ".Sun.bw")
        out_ilp_bw = os.path.join(args.assembly_dir, "hg38", g + ".ILP.bw")
        out_sun_ilp_bw = os.path.join(args.assembly_dir, "hg38", g + ".SUN_ILP.bw")
        out_ilp_counts = os.path.join(args.assembly_dir, "hg38", g + ".KmerCounts.bw")

        hg19_sun_counts = [os.path.join(args.output, g, "tracks", "{}.{}.UnfilteredSunModel.hg19.bedGraph".format(g, x)) for x in ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D"]]
        out_hg19 = [os.path.join(args.assembly_dir, "hg19", "{}.{}.UnfilteredSunModel.hg19.bedGraph".format(g, x)) for x in ["Notch2NL-A", "Notch2NL-B", "Notch2NL-C", "Notch2NL-D"]]

        if os.path.exists(bg) and os.path.exists(ilp) and os.path.exists(sun):
            system("bedGraphToBigWig {} {} {}".format(bg, args.hg38_chrom_sizes, out_sun_bw))
            system("wigToBigWig {} {} {}".format(ilp, args.hg38_chrom_sizes, out_ilp_bw))
            system("wigToBigWig {} {} {}".format(sun, args.hg38_chrom_sizes, out_sun_ilp_bw))
            system("wigToBigWig {} {} {}".format(ilp_counts, args.hg38_chrom_sizes, out_ilp_counts))
            paths[g] = [out_sun_bw, out_ilp_bw, out_sun_ilp_bw, out_ilp_counts]
        elif os.path.exists(bg) and os.path.exists(ilp):
            system("bedGraphToBigWig {} {} {}".format(bg, args.hg38_chrom_sizes, out_sun_bw))
            
            system("wigToBigWig {} {} {}".format(ilp, args.hg38_chrom_sizes, out_ilp_bw))
            system("wigToBigWig {} {} {}".format(ilp_counts, args.hg38_chrom_sizes, out_ilp_counts))
            paths[g] = [out_sun_bw, out_ilp_bw, None, out_ilp_counts]            
        elif os.path.exists(bg) and os.path.exists(sun):
            #probably amplicon
            system("bedGraphToBigWig {} {} {}".format(bg, args.hg38_chrom_sizes, out_sun_bw))
            system("wigToBigWig {} {} {}".format(sun, args.hg38_chrom_sizes, out_sun_ilp_bw))
            paths[g] = [out_sun_bw, None, out_sun_ilp_bw, None]            
        else:
            sys.stderr.write("Error: {} lacks tracks. Skipping.\n".format(g))
        for x,y in zip(hg19_sun_counts, out_hg19):
            system("bedGraphToBigWig {} {} {}".format(x, args.hg19_chrom_sizes, y))
        paths[g].append(out_hg19)
    buildTrackDb(args.assembly_dir, paths)

if __name__ == "__main__":
    main()
