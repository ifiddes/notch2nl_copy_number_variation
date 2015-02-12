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
    parser.add_argument("--chrom_sizes", help="hg38 chrom.sizes file", default="hg38.chrom.sizes")
    return parser.parse_args()


def startHub(d, name):
    with open(os.path.join(d, "hub.txt"), "w") as outf:
        outf.write("hub {}\nshortLabel {}\nlongLabel {}\ngenomesFile genomes.txt\nemail ian.t.fiddes@gmail.com\n".format(name, name + " Notch2NL", name + " Notch2NL"))
    with open(os.path.join(d, "genomes.txt"), "w") as outf:
        outf.write("genome hg38\ntrackDb hg38/trackDb.txt\ndefaultPos chr1:145987622-149723055\n")
    if not os.path.exists(os.path.join(d, "hg38")):
        os.mkdir(os.path.join(d, "hg38"))
    with open(os.path.join(d, "hg38", "Notch2NL.html"), "w") as outf:
         outf.write("Notch2NL {}\n".format(name))

def buildTrackDb(d, paths):
    with open(os.path.join(d, "hg38", "trackDb.txt"), "w") as outf:
        for g, [sun, ilp, sun_ilp] in paths.iteritems():
            outf.write("track {0}\ncontainer multiWig\nshortLabel {0}\nlongLabel {0}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\n\n".format(g))
            outf.write("track {0}_SUN\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}\n\n".format(g, os.path.basename(sun)))
            outf.write("track {0}_ILP\ncolor 242,148,176\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}\n\n\n".format(g, os.path.basename(ilp)))
            outf.write("track {0}_SUN_ILP\ncolor 209,203,205\nshortLabel {0}\nlongLabel {0}\nbigDataUrl {1}\ntype bigWig 0 4\nautoScale off\nvisibility full\nalwaysZero on\nyLineMark 2\nviewLimits 0:4\nyLineOnOff on\nmaxHeightPixels 100:75:50\nparent {0}\n\n\n".format(g, os.path.basename(sun_ilp)))


def main():
    args = parse_args()
    startHub(args.assembly_dir, args.name)

    genomes = [x for x in os.listdir(args.output)]
    paths = {}
    for g in genomes:
        bg = os.path.join(args.output, g, "tracks", g + ".UnfilteredSunModel.hg38.bedGraph")
        ilp = os.path.join(args.output, g, "tracks", g + ".ILP.wig")
        sun = os.path.join(args.output, g, "tracks", g + ".UnfilteredSunModel.hg38.SUN_ILP.wiggle")
        out_sun_bw = os.path.join(args.assembly_dir, "hg38", g + ".Sun.bw")
        out_ilp_bw = os.path.join(args.assembly_dir, "hg38", g + ".ILP.bw")
        out_sun_ilp_bw = os.path.join(args.assembly_dir, "hg38", g + ".SUN_ILP.bw")
        if os.path.exists(bg) and os.path.exists(ilp) and os.path.exists(sun):
            system("bedGraphToBigWig {} {} {}".format(bg, args.chrom_sizes, out_sun_bw))
            system("wigToBigWig {} {} {}".format(ilp, args.chrom_sizes, out_ilp_bw))
            system("wigToBigWig {} {} {}".format(sun, args.chrom_sizes, out_sun_ilp_bw))
            paths[g] = [out_sun_bw, out_ilp_bw, out_sun_ilp_bw]
        else:
            sys.stderr.write("Error: {} lacks tracks. Skipping.\n".format(g))
        

    buildTrackDb(args.assembly_dir, paths)

if __name__ == "__main__":
    main()
