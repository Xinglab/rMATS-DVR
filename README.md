# rMATS-DVR: rMATS calculation of Differential Variants of RNA

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
bam_calibration.py --bam test.bam --output /Path/to/output/test --genome hg19.fa --known dbSNP147.vcf
```	
	
###Required Parameters:
------------
	-h, --help       show this help message and exit

	--bam            Input bam file

	--output         Path and prefix of the output file

	--genome         genome sequence in fasta format

	--known          Known SNVs in vcf format



##2. Run rMATS-DVR to calculate Differential Variants of RNA
###Usage

```bash
rMATS-DVR.py [-h] [--sample1 SAMPLE1] [--sample2 SAMPLE2]
                    [--label LABEL] [--genome GENOME] [--known KNOWN]
                    [--output OUTPUT] [--minQ MINQ] [--minDP MINDP]
                    [--thread THREAD] [--diff DIFF]
```

























