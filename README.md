# rMATS-DVR: rMATS calculation of Differential Variants of RNA

##Requirements

1. Install Python 2.6.x or Python 2.7.x and corresponding versions of NumPy and SciPy. 
2. nstall Java.
3. Add the Python and Java directories to the $PATH environment variable.
4. Download Picard tool from https://broadinstitute.github.io/picard/.
5. Download GATK from https://software.broadinstitute.org/gatk/download/. (version 3.0~3.6 have been tested)

##Installation:

1. Unpack the downloaded tar ball. <br>
 - tar –zvxf rMATS-DVR-v1.0.tar.gz 
2. Create soft links of particular java programs from Picard and GATK into rMATS-DVR folder.<br>
 - cd rMATS-DVR-v1.0 <br>
 - ln -s  /path/to/picard/BuildBamIndex.jar ./ <br>
 - ln -s /path/to/picard/ReorderSam.jar ./ <br>
 - ln -s /path/to/picard/MarkDuplicates.jar ./ <br>
 - ln -s /path/to/picard/AddOrReplaceReadGroups.jar ./ <br>
 - ln -s /path/to/GATK/GenomeAnalysisTK.jar ./ <br>
3. Then the source code can be directly called from Python. <br>

###Required external files:
We highly recommend the users to use the genome sequence and dbSNP annotation from GTAK bundle, which can be downloaded from https://software.broadinstitute.org/gatk/download/bundle.

##1. Calibrate the bam files one by one.
###Usage
```bash
bam_calibration.py --bam sample.bam --output /Path/to/output/prefix --genome hg19.fa --known dbSNP147.vcf
```	
	
###Required Parameters:

	-h, --help       Show this help message and exit

	--bam            Input bam file

	--output         Path and prefix of the output file

	--genome         Genome sequence in fasta format

	--known          Known SNVs in vcf format



##2. Run rMATS-DVR to calculate Differential Variants of RNA
###Usage

```bash
rMATS-DVR.py --sample1 S1_rep1_calibrated.bam,S1_rep2_calibrated.bam,S1_rep3_calibrated.bam --sample2 S2_rep1_calibrated.bam,S2_rep2_calibrated.bam,S2_rep3_calibrated.bam --label S1,S2 --genome hg19.fa --known dbSNP147.vcf --output /Path/to/output/S1_vs_S2 [--minQ 20] [--minDP 5] [--thread 1] [--diff 0.0001]
```

###Required Parameters:

	-h, --help       Show this help message and exit

	--sample1        Bam files of sample 1, replicates are separated by comma
	
	--sample2        Bam files of sample 2, replicates are separated by comma
	
	--label          Lable of sample 1 and sample 2, separated by comma, e.g. Sample1,Sample2

	--output         Path and prefix of the output file

	--genome         Genome sequence in fasta format

	--known          Known SNVs in vcf format
	
###Optional Parameters:	
	--minQ MINQ      Minimum SNV quality [20]
	
	--minDP MINDP    Minimum mean read coverage of both samples [5]
	
	--thread         Number of processors [1]
	
	--diff           Required level difference between the two samples [0.0001]

