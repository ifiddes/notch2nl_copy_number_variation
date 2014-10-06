notch2nl_copy_number_variation
==============================
This pipeline is used to determine copy number of Notch2NL-A and Notch2NL-B from amplicon sequencing data.

Installation:
Clone this repo. In order to get the pipeline to work, you will need to copy a fasta of the hg19 assembly to the folder (it is too large for github). I have included the fasta index for this file as well as the chromsizes file used by kent tools.
You will want to name the fasta file hg19_ucsc.fa.

Usage:
The main script to run things is `amplicon_pipeline_wrapper_script.sh`. This script has one argument, a folder, which should contain sequencing results in the form of <id>_R1.fastq and <id>_R2.fastq. The important part is that the words R1 and R2 are present for forward and reverse reads.

This will then call `amplicon_pipeline.py` repeatedly for each sample in that folder.

`amplicon_pipeline.py` can be run individually as well, and has options to change the # of cores used for alignment as well the whitelist file and index. You may also change the --align flag to suppress realignments (generally for debugging purposes).

A full run of the pipeline will create an output folder in the remapped directory with the sample name. This folder will contain the raw bam file, the histogram based view, the ILP based view, and a bedGraph file that can be loaded into the genome browser (and whose coordinates are based on the CHM1 assembly).

If you want to include Notch2NL-C in the ILP plots, you will want to switch out the python scripts in `amplicon_pipeline.py` to use `integer_linear_programming_notch2nl_including_C.py`
