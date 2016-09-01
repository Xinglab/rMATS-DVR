# rMATS-DVR

Requirements:
1) Install Python 2.6.x or Python 2.7.x and corresponding versions of NumPy and SciPy.
 
2) Install Java. 

3) Add the Python and Java directories to the $PATH environment variable.

4) Download Picard tool from https://broadinstitute.github.io/picard/.

5) Download GATK from https://software.broadinstitute.org/gatk/download/. (version 3.0~3.6 have been tested)


Installation:


1) Unpack the downloaded tar ball.

	tar –zvxf rMATS-DVR-v1.0.tar.gz

2) Create soft links of particular java programs from Picard and GATK into rMATS-DVR folder. 

	cd rMATS-DVR-v1.0

	ln -s  /path/to/picard/BuildBamIndex.jar ./

	ln -s /path/to/picard/ReorderSam.jar ./

	ln -s /path/to/picard/MarkDuplicates.jar ./

	ln -s /path/to/picard/AddOrReplaceReadGroups.jar ./

	ln -s /path/to/GATK/GenomeAnalysisTK.jar ./


Usage:

1) Calibrate the bam files 


usage: bam_calibration.py [-h] [--bam BAM] [--output OUTPUT] [--genome GENOME]
                          [--known KNOWN]

Bam calibration for rMATS-DVR

optional arguments:

  -h, --help       show this help message and exit

  --bam BAM        Input bam file

  --output OUTPUT  Path and prefix of the output file

  --genome GENOME  genome sequence in fasta format

  --known KNOWN    Known SNVs in vcf format



