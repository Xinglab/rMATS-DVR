#!/usr/bin/env python

import argparse,os,sys,time,logging,datetime

startTime = time.time()

parser = argparse.ArgumentParser(description='Bam calibration for rMATS-DVR')
parser.add_argument('--bam',help='Input bam file')
parser.add_argument('--output',help='Path and prefix of the output file')
parser.add_argument('--genome',help='genome sequence in fasta format')
parser.add_argument('--known',help='Known SNVs in vcf format')

args = parser.parse_args()

bam=args.bam
label=args.output
genome=args.genome
known=args.known

directory='/'.join(sys.argv[0].split('/')[:-1])
if (directory!=''):
    directory+='/'

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=label+'.bamCalibration.log'+ str(datetime.datetime.now())+'.txt' ,
                    filemode='w')
  
# mark and index bam
com1='java -Xmx4g -jar '+directory+'ReorderSam.jar INPUT='+bam+' OUTPUT='+label+'_reordered.bam R='+genome
com2='java -Xmx4g -jar '+directory+'MarkDuplicates.jar INPUT='+label+'_reordered.bam OUTPUT='+label+'_dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE='+label+'_metrics.txt'
com3='java -Xmx4g -jar '+directory+'AddOrReplaceReadGroups.jar INPUT='+label+'_dedup.bam OUTPUT='+label+'_addrg.bam RGID='+label+' RGLB='+label+' RGPL=ILLUMINA RGPU=lane1 RGSM='+label
com4='java -Xmx4g -jar '+directory+'BuildBamIndex.jar INPUT='+label+'_addrg.bam'

# SplitNCigarReads
com5='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T SplitNCigarReads -R '+genome+' -I '+label+'_addrg.bam -o '+label+'_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'

# Base Quality Score Recalibration
com6='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T BaseRecalibrator -I '+label+'_split.bam -R '+genome+' -o '+label+'_recalibration_report.grp -knownSites '+known

com7='java -Xmx4g -jar '+directory+'GenomeAnalysisTK.jar -T PrintReads -R '+genome+' -I '+label+'_split.bam -BQSR '+label+'_recalibration_report.grp -o '+label+'_recalibration.bam -U ALLOW_N_CIGAR_READS'

logging.debug('Running command 1: '+com1+'\n')
os.system(com1)
logging.debug('Command 1 completed!\n')
logging.debug('Running command 2: '+com1+'\n')
os.system(com2)
logging.debug('Command 2 completed!\n')
logging.debug('Running command 3: '+com1+'\n')
os.system(com3)
logging.debug('Command 3 completed!\n')
logging.debug('Running command 4: '+com1+'\n')
os.system(com4)
logging.debug('Command 4 completed!\n')
logging.debug('Running command 5: '+com1+'\n')
os.system(com5)
logging.debug('Command 5 completed!\n')
logging.debug('Running command 6: '+com1+'\n')
os.system(com6)
logging.debug('Command 6 completed!\n')
logging.debug('Running command 7: '+com1+'\n')
os.system(com7)
logging.debug('Command 7 completed!\n')

logging.debug("Program ended")
currentTime = time.time()
runningTime = currentTime-startTime
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))

