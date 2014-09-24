#!/bin/bash


#produces an interleaved fastq from a bam in the region of notch2nl
BAM=$1
samtools view $BAM 1:120309986-145388461 -ub | samtools bamshuf -Ou /dev/stdin tmp | samtools bam2fq /dev/stdin > $BAM.fastq