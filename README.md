# rMATS-DVR: rMATS calculation of Differential variants of RNA

##Requirements
------------
1. Install Python 2.6.x or Python 2.7.x and corresponding versions of NumPy and SciPy. 
2. nstall Java.
3. Add the Python and Java directories to the $PATH environment variable.
4. Download Picard tool from https://broadinstitute.github.io/picard/.
5. Download GATK from https://software.broadinstitute.org/gatk/download/. (version 3.0~3.6 have been tested)

##Installation:
------------
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


##1. Calibrate the bam files one by one.
###Usage
```bash
	bam_calibration.py [--bam test.bam] [--output /Path/to/output/test] [--genome hg19.fa] [--known dbSNP147.vcf]
```bash


###Required Parameters:
------------
	  -h, --help       show this help message and exit

  	  --bam BAM        Input bam file

  --output OUTPUT  Path and prefix of the output file

  --genome GENOME  genome sequence in fasta format

  --known KNOWN    Known SNVs in vcf format

Optional Parameters:
------------	
	--o/--output:
		The output directory. The default is current directory
	--lib:
		The library type with choices of unstrand/first/second. The details are explained in the parameter of library-type in tophat2. The default is unstrand
	
	--read: 
		The sequencing strategy of producing reads with choices of P/S. The default is P
	--length: 
		The read length in nucleotide. The default length is 100
	--anchor: 
		The anchor length in nucleotide. The program will only count reads spanning junctions with at least this anchor length on each side. The default is 8
	--Cal: 
		Which  part of the program user choose to run, the choices are All/count/rMATS. All means run the whole program, count means only run the PI value calculation part, rMATS means only run the differential analysis of retained intron.  The default is All
	--RPKM: 
		A file providing the RPKM value for each sample, the first column is gene ID with the following column being the RPKM value for each sample. It is a required parameters to run the Density calculation
	--norm: 
		Total uniquely mapped reads for each library,each sample is seperated by comma, it is required to run the Density calculation
	--Clean: 
		true/false, whether to carry out PI_Density' calculation,The default is true
	--lim: 
		The minimum average number read per sample of the splice junction to be used in adjusting introns. The default value is 2
	--Comparison: 
		A file providing the sample pairs to calculate the differential RI level.The format should be column 1(name of comparions), column 2 (sample 1 order in the input file,replicates seperated by commas), column 3 (sample 2 order in the input file,replicates seperated by commas),column 4(type of PI value use to perform rMATS), column 5 (optional, if present as 'pool', the replicates are combined together). If absent, rMATS step will be skipped
	--analysis: 
		Type of rMATS analysis to perform. analysisType is either P or U. P is for paired analysis and U is for unpaired analysis. Default is U
	--c1: 
		The cutoff of splicing difference using Junction method. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001
	--c2: 
		The cutoff of splicing difference using Density method. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001
	--c3: 
		The cutoff of splicing difference using Density' method. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001
	--p: 
		The number of threads used to run rMATS. The default is 1;

Type of PI (Percent of Introns) Calculation:
------------	
	PI:
		Unspliced counts divided by transribed counts
	PI_Junction: 
		Inclusion counts divided by the sum of inclusion and  skipping junction counts
	PI_Density:
		The observed counts divided by the expected counts of the intron
	PI_Density': 
		The observed counts divided by the expected counts of the adjusted intron

Output list:
------------
	n = number of samples
	result:
		All of final result files are in result folder.
	counts_all_$type.txt store the inclusion and skipping counts for all of the samples
		column 1: Intron Id representing the chromosome position, start and end.
		column 2: Gene id
		column 3: Strand
		column 4: Chromosome name
		column 5: Start coordinate
		column 6: End coordinate
		column 7: Whether this intron was annotated in the gtf file as retained intron event.
		column 8: Whether this intron was overlapped with exon, or the 5' splice site was overlapped with exon or the 3' site was overlapped with exon or whether this intron is a simple intron
		column 9: Unspliced counts for all of the samples seperated by commas
		column 10: Transcribed counts for all of the samples seperated by commas
		column 11: Unspliced  length
		column 12: Transcribed length
		column 13: PI value for all of the samples seperated by commas
	rMATS_Result_$comparison_$type.txt store the differential RI level calculated by rMATS
		column 1: Intron Id representing the chromosome position, start and end.
		column 2: Gene id
		column 3: Strand
		column 4: Chromosome name
		column 5: Start coordinate
		column 6: End coordinate
		column 7: Whether this intron was annotated in the gtf file as retained intron event.
		column 8: Whether this intron was overlapped with exon, or the 5' splice site was overlapped with exon or the 3' site was overlapped with exon or whether this intron is a simple intron
		column 9: Unspliced counts for all replicates  of sample 1 seperated by commas
		column 10: Transcribed  counts for all replicates  of sample 1 seperated by commas
		column 11: Unspliced counts for all replicates  of sample 2 seperated by commas
		column 12: Transcribed counts for all replicates  of sample 2 seperated by commas
		column 13: Unspliced length
		column 14: Transcribed length
		column 15: p-value for differential PI level of the two samples
		column 16: FDR for differential RI level of the two samples
		column 17: PI level for sample1, replicates seperated by commas
		column 18: PI level for sample2, replicates seperated by commas
		column 19: The difference of PI level between sample1 and sample2, which is the result of average PI level of sample1 minus the average PI level of sampel2.

	test:
		A folder contains test files to run the program

	log.SQUID: Log file for running CARIE pipeline

	gtf_files:
		A folder contains different types of gtf files to run the program. Use mouse genome as examples.
	Mus_musculus.Ensembl.GRCm38.78.gtf: the ensemble gtf files. This file should be provided by user. 
	Exon_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains exons only
	Intron_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains intron only
	Intron_Annotated_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains the attributes whether the intron was annotated as retended introns in the original gtf files
	Intron_clean_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains the attributes whether the intron/5'Junction/3'Junction was overlapped with Exon and whether the intron is a simple intron. 

	counts:
		A folder contains all of the count files
	count_all.txt: a file contains the counts for all of the introns
		column 1:Intron Id representing the chromosome position, start and end.
		column 2:Gene id
		column 3:Strand
		column 4:Chromosome name
		column 5:Start coordinate
		column 6:End coordinate    
		column 7: Inclusion counts at 5' splice sites for sample 1
		column 8: Skipping counts at 5' splice sites for sample 1
		column 9: Inclusion counts at 3' splice sites for sample 1
		column 10: Skipping counts at 3' splice sites for sample 1
		column 11: Skipping counts of the intron for sample 1
		column 12: counts lying in the intron for sample 1
		column 13~6*(n+1): more counts for samples 2-n
	count_all_Density.txt: a file contains the observed counts and expected counts for all of the introns
		column 1: Intron id representing the chromosome position, start and end.
		column 2: The length of introns
		column 3~n+2: The observed counts
		column n+3~2n+2: The expected counts
	junction.txt: a file contains the spliced junction reads in the intron region
		column 1: Spliced junction reads id representing the chromosome position, start and end.
		column 2: The number of junction reads
		column 3~: The introns that the spliced junction reads reside
	Aintron.txt: a file contains the adjusted intron region
		column 1:Intron id representing the chromosome position, start and end
		column 2~: Adjusted intron region
	count_Clean_Density.txt: a file contains the count of adjusted intron
		column 1: Intron id representing the chromosome position, start and end
		column 2: Gene id
		column 3: The length of intron after adjusted intron
		column 4: Inclusion counts of at 5' splice sites of adjusted intron for sample 1
		column 5: Inclusion counts of at 3' splice sites of adjusted intron for sample 1
		column 6: Inclusion containing counts of adjusted intron for sample 1
		column 7~3(n+1): more counts for samples 2-n
	count_all_Clean_Density.txt: a file contains the observed counts and expected counts for all of the adjusted introns
		column 1: Intron id representing the chromosome position, start and end.
		column 2: The length of adjusted introns
		column 3~n+2: The observed counts of adjusted introns
		column n+3~2n+2: The expected counts of adjusted introns
		
	rMATS_files:
		A folder contains all of the rMATS input and output files
	rMATS_$comparison_$type.txt
    		The input file for running rMATS.
	rMATS_$comparison_$type folder
		The folder contains the result of rMATS output.
		
	cufflinks:
		A optional folder contains the result of cufflinks result and gene expression files for the squid run without gene expression file provided
	cufflinks_$n
		The cufflinks output folder of each folder
	gene_exp.txt
		The file is RPKM file that contains RPKM value for each gene.
		column 1: Gene ID
		column 2~n+1: RPKM value for samples 
FAQ
------------
1. Is the sample order in the Comparison file 1-based or 0-based?
A: The sample order in the Comparison file is 1-based.

