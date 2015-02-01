#!/usr/bin/env python
"""
Code by Adam Ewing.
"""
import argparse
import pysam
import vcf
import sys
import subprocess
import numpy as np
from os.path import basename
from re import sub

def is_somatic(rec):
    if str(rec.INFO.get('SS')).upper() in ['SOMATIC', '2']:
        return True

    if rec.INFO.get('SOMATIC'):
        if str(rec.INFO.get('SS')).upper() == 'LOH':
            return False
        return True

    if somatic_in_format(rec):
        return True

    return False

def somatic_in_format(rec):
    SS = []
    for sample in rec.samples:
        calldata = sample.data
        if 'SS' in calldata._fields:
            SS.append(calldata.SS)

    if '2' in SS or 2 in SS:
        return True
    return False

def basecount(bam,chrom,pos,prefixchr=False):
    baselist = [] 

    if 'chr' in chrom and not prefixchr:
        chrom = sub('chr','',chrom)
    if 'chr' not in chrom and prefixchr:
        chrom = 'chr' + chrom

    region = chrom + ":" + str(pos) + "-" + str(pos) 
    mpargs = ['samtools', 'mpileup', '-q','10','-r', region, bam.filename]

    mpileupstr=''
    p = subprocess.Popen(mpargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout.readlines():
        c = line.strip().split()
        if len(c) == 6:
           mpileupstr=c[4] 

    for base in list(mpileupstr):
        base = base.upper()
        if base in ['A', 'T', 'G', 'C']:
            baselist.append(base)

    bcount = {}
    for b in baselist:
        if b in bcount:
            bcount[b] += 1
        else:
            bcount[b] = 1

    return bcount

def main(args):
    invcf  = vcf.Reader(filename=args.vcffile)
    outvcf = vcf.Writer(sys.stdout, invcf)

    vtype = None
    if args.vtype is not None:
        assert args.vtype in ('SNV', 'INDEL', 'SV')
        vtype = args.vtype

    bam = pysam.Samfile(args.bamfile, 'rb')
    fa  = pysam.Fastafile(args.faidx)

    altfracs = []
    for rec in invcf:
        output = True
        altcount = 0
        bc = basecount(bam, rec.CHROM, rec.POS, prefixchr=args.chr)
        for alt in rec.ALT:
            if alt in bc.keys():
                altcount = int(bc[str(alt)])
                
            else:
                altcount = 0

        refbase = fa.fetch(rec.CHROM, int(rec.POS)-1, int(rec.POS)).upper()

        if refbase in rec.ALT:
            output = False

        if vtype == 'SNV' and (not rec.is_snp or (rec.is_snp and rec.INFO.get('VT') == 'LOH')):
            output = False

        if vtype == 'INDEL' and not rec.is_indel:
            output = False

        if vtype == 'SV' and not rec.is_sv:
            output = False

        if args.passonly and rec.FILTER:
            output = False

        if args.failonly and not rec.FILTER:
            output = False

        if args.somaticonly and not is_somatic(rec):
            output = False

        if args.germlineonly and is_somatic(rec):
            output = False

        totalcount = 0
        for b,c in bc.iteritems():
            totalcount += c
        if totalcount < int(args.minreads):
            output = False

        if output:
            altfrac = float(altcount)/float(totalcount) 
            altfracs.append(altfrac)
            rec.INFO['DP']  = totalcount
            rec.INFO['AF1'] = altfrac
            rec.INFO['AC1'] = altcount

            # add new info fields for clarity
            rec.INFO['ALTFRAC']  = altfrac 
            rec.INFO['ALTCOUNT'] = altcount
            rec.INFO['TOTCOUNT'] = totalcount
            rec.INFO['REFBASE']  = refbase
            outvcf.write_record(rec)
#        else:
#            if totalcount >= int(args.minreads):
#                altfracs.append(0.0)

#    ci = str(np.mean(altfracs)-2*np.std(altfracs)) + "-" + str(np.mean(altfracs)+2*np.std(altfracs))
#    sys.stderr.write("mean altfrac: " + str(np.mean(altfracs)) + " " + ci + " \n")
#    sys.stderr.write("median altfrac: " + str(np.median(altfracs)) + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Grab pileup columns from any number of BAM files corresponding to locations in a VCF file')
    parser.add_argument('-v', '--vcf', dest='vcffile', required=True, help='VCF file')
    parser.add_argument('-b', '--bam', dest='bamfile', required=True, help='BAM file')
    parser.add_argument('-i', '--index', dest='faidx', required=True, help='reference fasta (samtools indexed)')
    parser.add_argument('-m', '--minreads', dest='minreads', default=10, help='minimum reads (total coverage, default=10)')
    parser.add_argument('-c', '--context', dest='context', default=0, help='bases of context on either side of VCF entry')
    parser.add_argument('-t', '--vtype', dest='vtype', default=None, help='only include variants of vtype where vtype is SNV, INDEL, or SV')
    parser.add_argument('-p', '--passonly', action='store_true', default=False, help='only return PASS records')
    parser.add_argument('-f', '--failonly', action='store_true', default=False, help='only return non-PASS records')
    parser.add_argument('-s', '--somaticonly', action='store_true', default=False, help='only return somatic records')
    parser.add_argument('-g', '--germlineonly', action='store_true', default=False, help='only return germline records')
    parser.add_argument('--chr', action='store_true', default=False, help='add chr prefix to chromosome names')
    args = parser.parse_args()
    main(args)
