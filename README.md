# rMATS-DVR: rMATS calculation of Differential Variants of RNA

##Requirements

1. Install Python 2.6.x or Python 2.7.x and corresponding versions of NumPy and SciPy. 
2. Install Java.
3. Add the Python and Java directories to the $PATH environment variable.
4. Download Picard tool from https://broadinstitute.github.io/picard/.
5. Download GATK from https://software.broadinstitute.org/gatk/download/. (version 3.0~3.6 have been tested)
6. Download SAMtools (http://www.htslib.org/) and add it to $PATH environment variable.

##Installation:

1. Create soft links of particular java programs from Picard and GATK into rMATS-DVR folder.<br>
 - cd rMATS-DVR <br>
 - ln -s  /path/to/picard/BuildBamIndex.jar ./ <br>
 - ln -s /path/to/picard/ReorderSam.jar ./ <br>
 - ln -s /path/to/picard/MarkDuplicates.jar ./ <br>
 - ln -s /path/to/picard/AddOrReplaceReadGroups.jar ./ <br>
 - ln -s /path/to/GATK/GenomeAnalysisTK.jar ./ <br>
2. Then the source code can be directly called from Python. <br>

##Required and optional external files:

All the external files based on human hg19 genome can be downloaded from http://www.mimg.ucla.edu/faculty/xing/public_data/rMATS-DVR/hg19_resource.tar.gz

Alternatively, users can also prepare the external files under the following instructions:

1. Genome and Known SNV (required): we highly recommend the users to use the genome sequence and dbSNP annotation from GTAK bundle, which can be downloaded from https://software.broadinstitute.org/gatk/download/bundle. 
2. Known RNA editing sites: table delimited txt file with the first two columns are chromosome and coordinates. The other columns are ignored. Header is optional. Users can download the file from RADAR dababase (http://rnaedit.com/download/).
3. Genome-wide repeat elements: RepeatMasker Genomic Datasets downloaded from http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html. For example: hg19.fa.out.gz
4. Gene annotation: the gene annotaiton is in the UCSC format. We recommend users to download from UCSC. (http://hgdownload.soe.ucsc.edu/downloads.html#human). For example, one can download hg19 RefSeq gene from: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz


##Run rMATS-DVR in one step.

###Usage
In one step mode, rMATS-DVR will first calibrate bam files one by one and then calculate Differential Variants of RNA using all samples.

```bash
python rMATS-DVR.py --sample1 S1_rep_1.bam[,S1_rep_2.bam][,...,S1_rep_n.bam] --sample2 S2_rep_1.bam[,S2_rep_2.bam][,...,S2_rep_n.bam] --label S1,S2 --genome hg19.fa --known dbSNP147.vcf --output /Path/to/output/S1_vs_S2 [--editing RADAR2.txt] [--repeat repeats.txt] [--gene RefSeq.txt] [--minQ 20] [--minDP 5] [--thread 1] [--diff 0.0001] [--merge] [--ReadStranded] [--ReadPaired] [--skipBamCalibration] [--KeepTemp]
```

###Required Parameters:

	-h, --help              Show this help message and exit

	--sample1   <str>       Bam (or sam) files of sample 1, replicates are separated by comma
	
	--sample2   <str>       Bam (or sam) files of sample 2, replicates are separated by comma
	
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
	
	--ReadStranded          RNA-seq reads are Illumina strand-specific reads. Disable by default.
	
	--ReadPaired            RNA-seq reads are paired-end reads. Disable by default.
	
	--merge                 Merge the counts of all replicates. Enable by default when there are less than 2 replicates in at least one sample groups.
	
	--skipBamCalibration    Skip the step of calibrating bam files. Enable it when the input bam files have already been calibrated using bam_calibration.py (see below). Disable by default. 
	
	--KeepTemp		Keep the temporary files. Disable by default.
	

##Run rMATS-DVR in two steps.
When there are a large number of replicates, one step mode, which calibrate bam files one by one,  may take long time. In these cases, we recomend users to run bam calibration for all bam files parallely at the first step. Then the users can run rMATS-DVR.py with --skipBamCalibration.

Bam calibration
```bash
python bam_calibration.py --bam sample.bam --output /Path/to/output/prefix --genome hg19.fa --known dbSNP147.vcf
```	

###Required Parameters:

	-h, --help              Show this help message and exit

	--bam       <str>       Input bam (or sam) file

	--output    <str>       Path and prefix of the output file

	--genome    <str>       Genome sequence in fasta format

	--known     <str>       Known SNVs in vcf format


## Output

The final output files are in "Prefix_rMATS-DVR_results" folder, including "rMATS-DVR_Result.txt" and "rMATS-DVR_Result_summary.txt".
"rMATS-DVR_Result.txt" provides the SNV information, read counts, P value, FDR, gene location, and multiple annotations based on known databases. "rMATS-DVR_Result_summary.txt" summarizes the frequencies of all types of total SNVs and DVRs respectively. All other files are temporary files.

###1. rMATS-DVR_Result.txt
	 Ref_allele: reference allele.
	 Alt_allele: alternative allele.
	 RNA-seqStrand: RNA strand from which the RNA-seq reads are originated. Only valid when --ReadStranded is applied in rMATS-DVR.
	 Sample1_Alt: read counts of alternative allele in sample 1, replicates are separated by comma.
	 Sample1_Ref: read counts of reference allele in sample 1, replicates are separated by comma.
	 Sample2_Alt: read counts of alternative allele in sample 2, replicates are separated by comma.
	 Sample2_Ref: read counts of reference allele in sample 2, replicates are separated by comma.
	 Sample1_Alt_allele_fraction: fraction of alternative allele counts in sample 1, replicates are separated by comma.
	 Sample2_Alt_allele_fraction: fraction of alternative allele counts in sample 2, replicates are separated by comma.
	 Alt_allele_fraction_diff: average (Sample1_Alt_allele_fraction) - average (Sample2_Alt_allele_fraction).
	 Ref_onSense: reference allele on sense strand.
	 Alt_onSense: alternative allele on sense strand.
	 Location: location of the SNV on gene.
	 KnownSNV: rs ID of known SNV hit.
	 KnownRNAediting: boolean variable to show whether the SNV has a hit in known RNA editing database.
	 RepeatName: name of repeat element which covers the SNV.
	 RepeatName: family of repeat element which covers the SNV.

###2. rMATS-DVR_Result_summary.txt
	 Type (Ref-Alt) on sense strand: type of SNP in the format of reference allele-alternative allele on sense strand.
	 All SNV: frequency of each type of all called SNVs. 
	 DVR (FDR<0.05): frequency of each type of all SNVs with FDR <0.05. 

##Contacts

Yi Xing yxing@ucla.edu

Jinkai Wang jinkwang@ucla.edu

If you found a bug or mistake in this project, we would like to know about it. Before you send us the bug report though, please check the following:

Are you using the latest version? The bug you found may already have been fixed.
Check that your input is in the correct format and you have selected the correct options.
Please reduce your input to the smallest possible size that still produces the bug; we will need your input data to reproduce the problem, and the smaller you can make it, the easier it will be.
##Copyright and License Information

Copyright (C) 2016 University of California, Los Angeles (UCLA) Jinkai Wang and Yi Xing

Authors: Jinkai Wang and Yi Xing

This program is licensed with commercial restriction use license. Please see the attached LICENSE file for details.

