#!/usr/bin/env python

import argparse,os,sys,time,logging,datetime,re

startTime = time.time()

parser = argparse.ArgumentParser(description='rMATS-DVR (v1.0)')
parser.add_argument('--sample1',help='Bam files of sample 1, replicates are separated by comma')
parser.add_argument('--sample2',help='Bam files of sample 2, replicates are separated by comma')
parser.add_argument('--label',help='Lable of sample 1 and sample 2, separated by comma, e.g. Sample1,Sample2')
parser.add_argument('--genome',help='Genome sequence in fasta format')
parser.add_argument('--known',help='Known SNVs in vcf format')
parser.add_argument('--repeat', default='NA', help='Repeat elements annotation')
parser.add_argument('--editing', default='NA', help='Known RNA editing sites')
parser.add_argument('--gene', default='NA', help='Gene annotatione')
parser.add_argument('--output',help='Path and prefix of output file')
parser.add_argument('--minQ', default='20', help='Minimum SNV quality [20]')
parser.add_argument('--minDP', default='5', help='Minimum mean read coverage of both samples [5]')
parser.add_argument('--thread',  default='1', help='Number of processors [1]')
parser.add_argument('--diff',  default='0.0001', help='Required level difference between the two samples [0.0001]')
parser.add_argument('--merge', action='store_true', help='Merge the counts of all replicates. Enable by default when there are less than 2 replicates in at least one sample groups.')
parser.add_argument('--skipBamCalibration', action='store_true', help='Skip the step of calibrating bam files. Disable by default.')
parser.add_argument('--ReadStranded', action='store_true', help='RNA-seq reads are Illumina strand-specific reads. Disable by default.')
parser.add_argument('--ReadPaired', action='store_true', help='RNA-seq reads are paired-end reads. Disable by default.')
parser.add_argument('--KeepTemp', action='store_true', help='Keep tempory files. Disable by default.')

args = parser.parse_args()

sample1=args.sample1
sample2=args.sample2
genome=args.genome
known=args.known
repeatmask=args.repeat
knownediting=args.editing
geneanno=args.gene
output=args.output
minQ=args.minQ
minDP=args.minDP
thread=args.thread
diff=args.diff
lab1,lab2=args.label.split(',')
merge=args.merge
skip=args.skipBamCalibration
stranded=args.ReadStranded
paired=args.ReadPaired
keep=args.KeepTemp

ns1=len(sample1.split(','))
ns2=len(sample2.split(','))

if (ns1<2 or ns2<2):
    merge=True

directory='/'.join(sys.argv[0].split('/')[:-1])
if (directory!=''):
    directory+='/'
    
outdirectory='/'.join(output.split('/')[:-1])
if (outdirectory!=''):
    outdirectory+='/'

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=output+'_rMATS-DVR.log'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w')

samples=sample1.split(',')+sample2.split(',')
if (skip):
    allsample=' -I '.join(samples)
    newsamples=samples
else:
    allsample=[]
    newsamples=[]
    
if (not skip):
    logging.debug('Calibrating bam files\n')
    os.system('mkdir -p '+output+'_bam_calibration')
    for bam in samples:
        outbam=bam.split('/')[-1]
        logging.debug('Start calibrating '+outbam+'\n')
        calbam='.'.join(outbam.split('.')[:-1])
        allsample.append(output+'_bam_calibration/'+calbam+'_recalibration.bam')
        newsamples.append(output+'_bam_calibration/'+calbam+'_recalibration.bam')
        if (keep):
            com='python '+directory+'bam_calibration.py --bam '+bam+' --output '+output+'_bam_calibration/'+calbam+' --genome '+genome+' --known '+known+' --KeepTemp'
        else:
            com='python '+directory+'bam_calibration.py --bam '+bam+' --output '+output+'_bam_calibration/'+calbam+' --genome '+genome+' --known '+known
        os.system(com)
    allsample=' -I '.join(allsample)

com1='java -jar '+directory+'GenomeAnalysisTK.jar -T UnifiedGenotyper -R '+genome+' -I '+allsample+'  --dbsnp '+known+' -o '+output+'.vcf -stand_call_conf 0 -stand_emit_conf 0 --genotyping_mode DISCOVERY -nt '+thread
logging.debug('Running command 1: '+com1+'\n')
os.system(com1)

allbam_sep=newsamples
if (stranded and paired):
    logging.debug('Separate two end of bam files and pileup\n')
    allbam_sep=[]
    os.system('mkdir -p '+output+'_bam_seperate\n')
    for bam in newsamples:
        outbam=bam.split('/')[-1]
        bamname='.'.join(outbam.split('.')[:-1])
        os.system('samtools view -f 64 -b '+bam+' > '+ output+'_bam_seperate/'+bamname+'_firstEnd.bam' )
        os.system('samtools view -f 128 -b '+bam+' > '+ output+'_bam_seperate/'+bamname+'_secondEnd.bam' )
        allbam_sep.append(output+'_bam_seperate/'+bamname+'_firstEnd.bam')
        allbam_sep.append(output+'_bam_seperate/'+bamname+'_secondEnd.bam')
os.system('samtools mpileup -B -d 100000 -f '+genome+' -l '+output+'.vcf -q 30 -Q 17 -o '+output+'.pileup '+' '.join(allbam_sep))
if (not keep):
    os.system('rm -rf '+output+'_bam_calibration')
    os.system('rm -rf '+output+'_bam_seperate')

if (stranded and paired):
    com2='python '+directory+'vcf_to_mats_input.py '+output+'.vcf '+output+'.inc.txt '+sample1+' '+sample2+' '+minQ+' '+minDP+' T '+output+'.pileup T'
elif(stranded and not paired):
    com2='python '+directory+'vcf_to_mats_input.py '+output+'.vcf '+output+'.inc.txt '+sample1+' '+sample2+' '+minQ+' '+minDP+' F '+output+'.pileup T'
else:
    com2='python '+directory+'vcf_to_mats_input.py '+output+'.vcf '+output+'.inc.txt '+sample1+' '+sample2+' '+minQ+' '+minDP+' F '+output+'.pileup F'
if (not merge):
    com3='python '+directory+'rMATS_unpaired.py '+output+'.inc.txt '+output+'_rMATS-DVR_results '+thread+' '+diff
else:
    com3='python '+directory+'MATS_LRT.py '+output+'.inc.txt '+output+'_rMATS-DVR_results '+thread+' '+diff
com4='python '+directory+'FDR.py '+output+'_rMATS-DVR_results/rMATS_Result_P.txt '+output+'_rMATS-DVR_results/rMATS_Result_FDR.txt'
com5='python '+directory+'snv_annotation.py --input '+output+'_rMATS-DVR_results/rMATS_Result_FDR.txt --output '+output+'_rMATS-DVR_results/rMATS-DVR_Result.txt --summary '+output+'_rMATS-DVR_results/rMATS-DVR_Result_summary.txt --label1 '+lab1+' --label2 '+lab2+' --snp '+known+' --repeat '+repeatmask+' --editing '+knownediting+' --gene '+geneanno


logging.debug('Running command 2: '+com2+'\n')
os.system(com2)
logging.debug('Running command 3: '+com3+'\n')
os.system('mkdir -p '+output+'_rMATS-DVR_results')
os.system(com3)
logging.debug('Running command 4: '+com4+'\n')
os.system(com4)
logging.debug('Running command 5: '+com5+'\n')
os.system(com5)

logging.debug("Program ended")
currentTime = time.time()
runningTime = currentTime-startTime
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))


