import sys

vcf_input=sys.argv[1]
output=open(sys.argv[2], 'w')
sample1=sys.argv[3].split(',')
sample2=sys.argv[4].split(',')
minQ=int(sys.argv[5])
minDP=int(sys.argv[6])

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
        site=':'.join([chrom, pos, ref, alt, qual, '.'])
        ctr_refdp=[]
        ctr_altdp=[]
        hyp_refdp=[]
        hyp_altdp=[]
        n_ctr=0
        n_hyp=0
        for i in range(9, 9+ns1):
            b=a[i].split(':')
            if (b[0]=='./.'):
                ref_dp='0'
                alt_dp='0'
                n_ctr+=1
            else:  
                ref_dp, alt_dp=b[1].split(',')
            ctr_refdp.append(ref_dp)
            ctr_altdp.append(alt_dp)
        for i in range(9+ns1,9+ns1+ns2):
            b=a[i].split(':')
            if (b[0]=='./.'):
                ref_dp='0'
                alt_dp='0'
                n_hyp+=1
            else:
                ref_dp, alt_dp=b[1].split(',')
            hyp_refdp.append(ref_dp)
            hyp_altdp.append(alt_dp)
        if (n_ctr>1 or n_hyp>1):
            continue
        ctr_dp=[]
        hyp_dp=[]
        for i in range(ns1):
            ctr_dp.append(int(ctr_refdp[i])+int(ctr_altdp[i]))
        for i in range(ns2):
            hyp_dp.append(int(hyp_refdp[i])+int(hyp_altdp[i]))
        if (sum(ctr_dp)< ns1*minDP or sum(hyp_dp)<ns2*minDP):
            continue
        output.write('\t'.join([site, ','.join(ctr_altdp), ','.join(ctr_refdp), ','.join(hyp_altdp), ','.join(hyp_refdp), '100','100'])+'\n')

        

