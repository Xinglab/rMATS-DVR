# rMATS-DVR: rMATS calculation of Differential Variants of RNA

##Requirements

1. Install Python 2.6.x or Python 2.7.x and corresponding versions of NumPy and SciPy. 
2. Install Java.
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

##Required or optional external files:

All the external files of human hg19 genome can be downloaded from http://www.mimg.ucla.edu/faculty/xing/public_data/rMATS-DVR/hg19_resource.tar.gz

Alternatively, users can also prepare the external files under the following instructions:

1. Genome and Known SNV (required): we highly recommend the users to use the genome sequence and dbSNP annotation from GTAK bundle, which can be downloaded from https://software.broadinstitute.org/gatk/download/bundle. 
2. Known RNA editing sites: table delimited txt file with the first two columns are chromosome and coordinates. The other columns are ignored. Header is optional. Users can download the file from RADAR dababase (http://rnaedit.com/download/).
3. Genome-wide repeat elements: RepeatMasker Genomic Datasets downloaded from http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html. For example: hg19.fa.out.gz
4. Gene annotation: the gene annotaiton is in the UCSC format. We recommend users to download from UCSC. (http://hgdownload.soe.ucsc.edu/downloads.html#human). For example, one can download hg19 RefSeq gene from: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz


##Run rMATS-DVR in one step.

###Usage
In one step mode, rMATS-DVR will first calibrate bam files one by one and then calculate Differential Variants of RNA using all samples.

```bash
rMATS-DVR.py --sample1 S1_rep_1.bam[,S1_rep_2.bam][,...,S1_rep_n.bam] --sample2 S2_rep_1.bam[,S2_rep_2.bam][,...,S2_rep_n.bam] --label S1,S2 --genome hg19.fa --known dbSNP147.vcf --output /Path/to/output/S1_vs_S2 [--editing RADAR2.txt] [--repeat repeats.txt] [--gene RefSeq.txt] [--minQ 20] [--minDP 5] [--thread 1] [--diff 0.0001] [--merge] [--skipBamCalibration]
```

###Required Parameters:

	-h, --help              Show this help message and exit

	--sample1   <str>       Bam files of sample 1, replicates are separated by comma
	
	--sample2   <str>       Bam files of sample 2, replicates are separated by comma
	
	--label     <str>       Lable of sample 1 and sample 2, separated by comma, e.g. Sample1,Sample2

	--output    <str>       Path and prefix of the output file

	--genome    <str>       Genome sequence in fasta format

	--known     <str>       Known SNVs in vcf format
	
###Optional Parameters:

	--editing   <str>       Known RNA editing sites
	
	--repeat    <str>       Repeat elements annotation
	
	--gene      <str>       Gene annotation

	--minQ      <int>       Minimum SNV quality [20]
	
	--minDP     <int>       Minimum mean read coverage of both samples [5]
	
	--thread    <int>       Number of processors [1]
	
	--diff      <float>     Required level difference between the two samples [0.0001]
	
	--merge                 Merge the counts of all replicates. Enable by default when there are less than 2 replicates in at least one sample groups.
	
	--skipBamCalibration    Skip the step of calibrating bam files. Enable it when the input bam files have already been calibrated using bam_calibration.py (see below). Disable by default. 



##Run rMATS-DVR in two steps.
When there are a large number of replicates, one step mode, which calibrate bam files one by one,  may takes long time. In these cases, we recomend users to run bam calibration for all bam files parallely at the first step. Then the users can run rMATS-DVR.py with --skipBamCalibration.

Bam calibration
```bash
bam_calibration.py --bam sample.bam --output /Path/to/output/prefix --genome hg19.fa --known dbSNP147.vcf
```	

###Required Parameters:

	-h, --help              Show this help message and exit

	--bam       <str>       Input bam file

	--output    <str>       Path and prefix of the output file

	--genome    <str>       Genome sequence in fasta format

	--known     <str>       Known SNVs in vcf format

