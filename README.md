# rMATS-DVR 

Requirements: <br>
1) Install Python 2.6.x or Python 2.7.x and corresponding versions of NumPy and SciPy.<br>
2) Install Java. <br>
3) Add the Python and Java directories to the $PATH environment variable. <br>
4) Download Picard tool from https://broadinstitute.github.io/picard/. <br>
5) Download GATK from https://software.broadinstitute.org/gatk/download/. (version 3.0~3.6 have been tested) <br>
<br>
<br>
Installation:<br>
1) Unpack the downloaded tar ball. <br>
	tar –zvxf rMATS-DVR-v1.0.tar.gz <br>
2) Create soft links of particular java programs from Picard and GATK into rMATS-DVR folder. <br>
	cd rMATS-DVR-v1.0 <br>
	ln -s  /path/to/picard/BuildBamIndex.jar ./ <br>
	ln -s /path/to/picard/ReorderSam.jar ./ <br>
	ln -s /path/to/picard/MarkDuplicates.jar ./ <br>
	ln -s /path/to/picard/AddOrReplaceReadGroups.jar ./ <br>
	ln -s /path/to/GATK/GenomeAnalysisTK.jar ./ <br>
<br>
<br>
Usage: <br>
1) Calibrate the bam files <br>
usage: bam_calibration.py [-h] [--bam BAM] [--output OUTPUT] [--genome GENOME] [--known KNOWN] <br>
Bam calibration for rMATS-DVR <br>
<br>
optional arguments: <br>
  -h, --help       show this help message and exit <br>
  --bam BAM        Input bam file <br>
  --output OUTPUT  Path and prefix of the output file <br>
  --genome GENOME  genome sequence in fasta format <br>
  --known KNOWN    Known SNVs in vcf format <br>


