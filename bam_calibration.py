#!/usr/bin/env python

import argparse,os,sys,time,logging,datetime

startTime = time.time()

parser = argparse.ArgumentParser(description='Bam calibration for rMATS-DVR')
parser.add_argument('--bam',help='Input bam file')
parser.add_argument('--output',help='Path and prefix of the output file')
parser.add_argument('--genome',help='genome sequence in fasta format')
parser.add_argument('--known',help='Known SNVs in vcf format')
parser.add_argument('--KeepTemp', action='store_true', help='Keep tempory files. Disable by default.')

args = parser.parse_args()

bam=args.bam
label=args.output
genome=args.genome
known=args.known
keep=args.KeepTemp

directory='/'.join(sys.argv[0].split('/')[:-1])
if (directory!=''):
    directory+='/'

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=label+'.bamCalibration.log'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w')
  
# mark and index bam
com1='java -Xmx4g -jar '+directory+'picard.jar ReorderSam INPUT='+bam+' OUTPUT='+label+'_reordered.bam S=true R='+genome
com2='java -Xmx4g -jar '+directory+'picard.jar AddOrReplaceReadGroups INPUT='+label+'_reordered.bam OUTPUT='+label+'_addrg.bam RGID='+label+' RGLB='+label+' RGPL=ILLUMINA RGPU=lane1 RGSM='+label
com3='java -Xmx4g -jar '+directory+'picard.jar MarkDuplicates INPUT='+label+'_addrg.bam OUTPUT='+label+'_dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE='+label+'_metrics.txt'

# SplitNCigarReads
com5='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T SplitNCigarReads -R '+genome+' -I '+label+'_dedup.bam -o '+label+'_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'

# Base Quality Score Recalibration
com6='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T BaseRecalibrator -I '+label+'_split.bam -R '+genome+' -o '+label+'_recalibration_report.grp -knownSites '+known

com7='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T PrintReads -R '+genome+' -I '+label+'_split.bam -BQSR '+label+'_recalibration_report.grp -o '+label+'_recalibration.bam -U ALLOW_N_CIGAR_READS'

logging.debug('Running command 1: '+com1+'\n')
os.system(com1)
logging.debug('Running command 2: '+com2+'\n')
os.system(com2)
if (not keep):
    os.system('rm -f '+label+'_reordered.bam')
logging.debug('Running command 3: '+com3+'\n')
os.system(com3)
if (not keep):
    os.system('rm -f '+label+'_addrg.bam')
logging.debug('Running command 5: '+com5+'\n')

os.system(com5)
if (not keep):
    os.system('rm -f '+label+'_dedup.bam')
logging.debug('Running command 6: '+com6+'\n')
os.system(com6)
logging.debug('Running command 7: '+com7+'\n')
os.system(com7)
if (not keep):
    os.system('rm -f '+label+'_split.bam')
    os.system('rm -f '+label+'_recalibration_report.grp')
    os.system('rm -f '+label+'_split.bai')
    os.system('rm -f '+label+'_metrics.txt')
    
logging.debug("Program ended")
currentTime = time.time()
runningTime = currentTime-startTime
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))

