#!/bin/bash

fwd=$1
rev=$2
out=$3
cores=$4
index=$5

bwa mem -t $cores -v 1 $index $fwd $rev | samtools view -@ $cores -bS - | samtools sort -@ $cores - $out.sorted

python filter_and_fix_coords.py $out.sorted.bam $out.remapped.sorted.bam

samtools index $out.remapped.sorted.bam

rm $out.sorted.bam

#samtools mpileup -uf $index $out.remapped.sorted.bam | bcftools view -bvcg - | bcftools view - > $out.vcf

vcf_n=vcfs/N.vcf.gz
vcf_a=vcfs/A.vcf.gz
vcf_b=vcfs/B.vcf.gz
vcf_c=vcfs/C.vcf.gz
vcf_d=vcfs/D.vcf.gz
vcf_ab=vcfs/AB.vcf.gz

python getsitecoverage.py -t SNV -b $out.remapped.sorted.bam -v $vcf_n -i hg19_ucsc.fa --chr > $out.N.vcf
python getsitecoverage.py -t SNV -b $out.remapped.sorted.bam -v $vcf_a -i hg19_ucsc.fa --chr > $out.A.vcf
python getsitecoverage.py -t SNV -b $out.remapped.sorted.bam -v $vcf_b -i hg19_ucsc.fa --chr > $out.B.vcf
python getsitecoverage.py -t SNV -b $out.remapped.sorted.bam -v $vcf_ab -i hg19_ucsc.fa --chr > $out.AB.vcf
python getsitecoverage.py -t SNV -b $out.remapped.sorted.bam -v $vcf_c -i hg19_ucsc.fa --chr > $out.C.vcf
python getsitecoverage.py -t SNV -b $out.remapped.sorted.bam -v $vcf_d -i hg19_ucsc.fa --chr > $out.D.vcf
