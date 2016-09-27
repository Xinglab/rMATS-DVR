import sys,re

vcf_input=sys.argv[1]
output=open(sys.argv[2], 'w')
sample1=sys.argv[3].split(',')
sample2=sys.argv[4].split(',')
minQ=int(sys.argv[5])
minDP=int(sys.argv[6])
paired=sys.argv[7]   # T or F, whether the reads are paired, for unstranded, just treat it like single end (use F)
pileup=sys.argv[8]  # output of samtools mpileup
stranded=sys.argv[9]  # T or F, whether the reads are stranded.

if (stranded=='F'):
    paired='F'

def basecount (ref, seq, end):
    f_types = {'A':0, 'G':0,'C':0,'T':0}
    r_types= {'A':0, 'G':0,'C':0,'T':0}
    ref=ref.upper()
    i=0
    while i <len(seq):
        if (seq[i]=='+' or seq[i]=='-'):
            num=re.match(r"(\d+)", seq[i+1:]).group(0)
            i+=len(num)+int(num)
        elif ((seq[i] in ['A', 'G', 'C', 'T'] and end==1) or (seq[i] in ['a', 'g', 'c', 't'] and end==2)):
            r_types[seq[i].upper()]+=1
        elif((seq[i] in ['A', 'G', 'C', 'T'] and end==2) or (seq[i] in ['a', 'g', 'c', 't'] and end==1)):
            f_types[seq[i].upper()]+=1
        elif ((seq[i]=='.' and end==1) or (seq[i]==',' and end==2)):
            r_types[ref]+=1
        elif ((seq[i]==',' and end==1) or (seq[i]=='.' and end==2)):
            f_types[ref]+=1
        else:
            pass
        i+=1
    return [f_types, r_types]


f_count={}
r_count={}
for line in open(pileup):
            line=line.rstrip('\n\r')
            a=line.split('\t')
            chrom=a[0]
            site=a[1]
            ref=a[2].upper()
            f_count[chrom+':'+site]=[]
            r_count[chrom+':'+site]=[]
            if (paired=='F'):
                for i in range(4, len(a), 3):
                    count=basecount(ref, a[i], 1)
                    f_count[chrom+':'+site].append(count[0])
                    r_count[chrom+':'+site].append(count[1])
            elif (paired=='T'):
                for i in range(4, len(a), 6):
                    first=basecount(ref, a[i],1)
                    second=basecount(ref, a[i+3], 2)
                    f_one={}
                    r_one={}
                    for base in ['A', 'G', 'C', 'T']:
                        f_one[base]=first[0][base]+second[0][base]
                        r_one[base]=first[1][base]+second[1][base]
                    f_count[chrom+':'+site].append(f_one)
                    r_count[chrom+':'+site].append(r_one)
all_count=[f_count, r_count]

output.write('\t'.join(['ID',	'IC_SAMPLE_1', 'SC_SAMPLE_1', 'IC_SAMPLE_2', 'SC_SAMPLE_2', 'IncFormLen', 'SkipFormLen'])+'\n')

ns1=len(sample1)
ns2=len(sample2)

for line in open(vcf_input):
        line=line.rstrip('\n\r')
        if (line[0]=='#'):
            continue
        a=line.split('\t')
        chrom=a[0]
        pos=a[1]
        ref=a[3]
        alt=a[4]
        qual=a[5]
        if (float(qual)<minQ):
            continue
        if (len(ref)!=1 or len(alt)!=1):
            continue
        if (a[8]!='GT:AD:DP:GQ:PL'):
            continue
        if (stranded=='T'):
            strands=['+', '-']
            for m in range(2):
                site=':'.join([chrom, pos, ref, alt, qual, strands[m]])
                ctr_refdp=[]
                ctr_altdp=[]
                hyp_refdp=[]
                hyp_altdp=[]
                ctr_dp=[]
                hyp_dp=[]
                for i in range(ns1):
                        ctr_refdp.append(str(all_count[m][chrom+':'+pos][i][ref]))
                        ctr_altdp.append(str(all_count[m][chrom+':'+pos][i][alt]))
                for i in range(ns1, ns1+ns2):
                        hyp_refdp.append(str(all_count[m][chrom+':'+pos][i][ref]))
                        hyp_altdp.append(str(all_count[m][chrom+':'+pos][i][alt]))
                for i in range(ns1):
                    ctr_dp.append(int(ctr_refdp[i])+int(ctr_altdp[i]))
                for i in range(ns2):
                    hyp_dp.append(int(hyp_refdp[i])+int(hyp_altdp[i]))
                if (sum(ctr_dp)< ns1*minDP or sum(hyp_dp)<ns2*minDP):
                    continue
                output.write('\t'.join([site, ','.join(ctr_altdp), ','.join(ctr_refdp), ','.join(hyp_altdp), ','.join(hyp_refdp), '100','100'])+'\n')
        elif (stranded=='F'):
                site=':'.join([chrom, pos, ref, alt, qual, '.'])
                ctr_refdp=[]
                ctr_altdp=[]
                hyp_refdp=[]
                hyp_altdp=[]
                ctr_dp=[]
                hyp_dp=[]
                for i in range(ns1):
                        ctr_refdp.append(str(all_count[0][chrom+':'+pos][i][ref]+ all_count[1][chrom+':'+pos][i][ref]))
                        ctr_altdp.append(str(all_count[0][chrom+':'+pos][i][alt]+ all_count[1][chrom+':'+pos][i][alt]))
                for i in range(ns1, ns1+ns2):
                        hyp_refdp.append(str(all_count[0][chrom+':'+pos][i][ref]+ all_count[1][chrom+':'+pos][i][ref]))
                        hyp_altdp.append(str(all_count[0][chrom+':'+pos][i][alt]+ all_count[1][chrom+':'+pos][i][alt]))
                for i in range(ns1):
                    ctr_dp.append(int(ctr_refdp[i])+int(ctr_altdp[i]))
                for i in range(ns2):
                    hyp_dp.append(int(hyp_refdp[i])+int(hyp_altdp[i]))
                if (sum(ctr_dp)< ns1*minDP or sum(hyp_dp)<ns2*minDP):
                    continue
                output.write('\t'.join([site, ','.join(ctr_altdp), ','.join(ctr_refdp), ','.join(hyp_altdp), ','.join(hyp_refdp), '100','100'])+'\n')

        

