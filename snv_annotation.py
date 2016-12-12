import sys,argparse,re

"""
argv1:  input file
argv2: output file
argv3: Results summary table
argv4: Label of sample1
argv5: Label of sample2
argv6: dbSNP vcf
"""

parser = argparse.ArgumentParser(description='SNV annotation')
parser.add_argument('--input',help='Input file')
parser.add_argument('--output',help='Output file')
parser.add_argument('--summary',help='Results summary file')
parser.add_argument('--label1',help='Label of sample 1')
parser.add_argument('--label2',help='Label of sample 2')
parser.add_argument('--snp',help='Known SNVs in vcf format')
parser.add_argument('--repeat', default='NA', help='Repeat annotation')
parser.add_argument('--editing', default='NA', help='Known RNA editing sites')
parser.add_argument('--gene', default='NA', help='Gene annotation in gtf format')
args = parser.parse_args()

repeatmask=args.repeat
knownediting=args.editing
geneanno=args.gene

lab1=args.label1
lab2=args.label2

complement={'A':'T', 'T':'A', 'C':'G', 'G':'C', '-':'-', ',':','}

def rev_comp (seq):
            seq=seq.upper()
            newseq=''
            for base in seq:
                        newseq+=complement[base]
            newseq=newseq[::-1]
            return newseq

def exon_cds (start, end, cdsstart,cdsend):
            for i in range(len(start)):
                if (start[i]<=cdsstart<end[i]):              
                    m=i
                elif(cdsstart==end[i]):
                    #print ['cdsstart==end[i]', start, end, cdsstart,cdsend]
                    m=i+1
                if (start[i]<=cdsend<=end[i]):
                    n=i
                    break
            #print['exon_cds', start, end, cdsstart,cdsend]
            changedstart=[cdsstart]+start[m+1:n+1]
            changedend=end[m:n]+[cdsend]
            utr5start=start[:m+1]
            utr5end=end[:m]+[cdsstart]
            utr3start=[cdsend]+start[n+1:]
            utr3end=end[n:]
            cdslength=0
            utr5length=0
            utr3length=0
            for i in range(len(changedstart)):
                cdslength+=changedend[i]-changedstart[i]
            for i in range(len(utr5start)):
                utr5length+=utr5end[i]-utr5start[i]
            for i in range(len(utr3start)):
                utr3length+=utr3end[i]-utr3start[i]
            return  [changedstart, changedend,cdslength, utr5start, utr5end, utr5length,utr3start, utr3end, utr3length]

txbins={}
gene={}        #key is tx
if geneanno!='NA':
      for line in open(geneanno):
            line=line.rstrip('\n\r')
            a=line.split('\t')
            if (line[0]=='#'):
                continue
            chrom=a[2]
            tx=a[1] 
            #genename=a[14]
            genename=a[12]
            strand=a[3]
            txstart=a[9].split(',')[:-1]
            txend=a[10].split(',')[:-1]
            for i in range(len(txstart)):
                        txstart[i]=int(txstart[i])
                        txend[i]=int(txend[i])
            cdsstart=int(a[6])
            cdsend=int(a[7])
            change=exon_cds(txstart, txend, cdsstart, cdsend)
            gene[tx]=[chrom]+change[:]+[cdsstart]+[cdsend]+[strand]+[genename]
            genestart=int(a[4])
            geneend=int(a[5])
            bin1=genestart/1000000
            bin2=(geneend-1)/1000000
            if (not txbins.has_key(chrom)):
                txbins[chrom]={}
            for onebin in range(bin1, bin2+1):
                if (not txbins[chrom].has_key(onebin)):
                    txbins[chrom][onebin]=[tx]
                else:
                    txbins[chrom][onebin].append(tx)

allsites={}
for line in open(args.input):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    if (a[0]=='ID'):
        continue
    chr, site, ref, alt, qual, strand=a[0].split(':')
    allsites[':'.join([chr, site])]=[ref, alt, qual, strand]

print 'All sites loaded'

dbsnp={}
for line in open(args.snp):
    line=line.rstrip('\n\r')
    if (line[0]=='#'):
        continue
    a=line.split('\t')
    chrom=a[0]
    site=a[1]
    if (allsites.has_key(chrom+':'+site)):
        rs=a[2]
        ref=a[3]
        alt=a[4]
        dbsnp[chrom+':'+site]=[rs, ref, alt]
    
print 'dbSNP loaded'

radar2={}
if (knownediting!='NA'):
            for line in open(knownediting):
                line=line.rstrip('\n\r')
                a=line.split('\t')
                if (a[0]=='chromosome'):
                    continue
                if (allsites.has_key(a[0]+':'+a[1])):
                    radar2[a[0]+':'+a[1]]=1

print 'Known RNA editing loaded'

repeat={}
if repeatmask!='NA':
        for line in open(repeatmask):
                line=line.rstrip('\n\r')
                a=re.split(r"\s+", line)
                a=a[1:]
                if (len(a)<5 or len(a[4])<3 or  a[4][:3]!='chr'):
                    continue
                chrom=a[4]
                start=int(a[5])-1
                end=int(a[6])
                repname=a[9]
                repfam=a[10]
                bin1=start/100000
                bin2=(end-1)/100000
                if (not repeat.has_key(chrom)):
                    repeat[chrom]={}
                for onebin in range(bin1, bin2+1):
                    if(not repeat[chrom].has_key(onebin)):
                        repeat[chrom][onebin]=[[start, end, repname, repfam]]  
                    else:
                        repeat[chrom][onebin].append([start, end, repname, repfam])

print 'repeats loaded'
 
output=open(args.output, 'w')

types=['A-C', 'A-G', 'A-T', 'C-A', 'C-G', 'C-T', 'G-A', 'G-C', 'G-T', 'T-A', 'T-C', 'T-G']
all_type={}
dvr_type={}
snp_type={}
editing_type={}
novel_type={}
for item in types:
    all_type[item]=0
    dvr_type[item]=0
    snp_type[item]=0
    editing_type[item]=0
    novel_type[item]=0

output.write('\t'.join(['Chrom', 'Site', 'Ref_allele', 'Alt_allele', 'RNA-seqStrand', 'SNV_quality', lab1+'_Alt', lab1+'_Ref', lab2+'_Alt', lab2+'_Ref', 'Pvalue','FDR', lab1+'_Alt_allele_fraction', lab2+'_Alt_allele_fraction', 'Alt_allele_fraction_diff',
            'Genename', 'strand', 'Ref_onSense', 'Alt_onSense', 'Location', 'KnownSNV', 'KnownRNAediting', 'RepeatName', 'RepeatFamily'])+'\n')
for line in open(args.input):
    line=line.rstrip('\n\r')
    a=line.split('\t')
    if (a[0]=='ID'):     
        continue
    chr, site, ref, alt, qual, RNAstrand=a[0].split(':')
    site=int(site)-1
    ref=ref.upper()
    alt=alt.upper()
    tx=[]
    loctype=[]
    genename=[]
    #ensg=[]
    strands=[]
    sitebin=site/1000000
    #print [chr, sitebin,  len(txbins[chr]), len(txbins[chr][sitebin])]
    if (txbins.has_key(chr) and txbins[chr].has_key(sitebin)):
        for key in txbins[chr][sitebin]:
            intron=0
            cdsstart=[]
            cdsend=[]
            length=0
            cdsseq=''
            location=''
            utr5start=[]
            utr5end=[]
            utr3start=[]
            utr3end=[]
            txchr=gene[key][0]
            txstrand=gene[key][12]
            txgenename=gene[key][13]
            if(gene[key][4][0]<=site<gene[key][8][-1] and (RNAstrand=='.'  or RNAstrand==txstrand)):
                if (not gene[key][13] in genename):
                    genename.append(gene[key][13])
                    strands.append(txstrand)
            if (gene[key][4][0]<=site<gene[key][5][-1] and (RNAstrand=='.'  or RNAstrand==txstrand)):
                utr5start=gene[key][4]
                utr5end=gene[key][5]
                tx.append(key)
                for j in range(len(utr5start)):
                    if (utr5start[j]<=site<utr5end[j]):
                        if (gene[key][1]==gene[key][2]):
                            location='ncRNA'
                        elif (txstrand=='+'):
                            location='5UTR'
                        elif (txstrand=='-'):
                            location='3UTR'
                        break
                    elif(site>=utr5end[j]):
                        pass
                    elif(site<utr5start[j]):
                        intron=1
                        location='Intron'
                        break
                loctype.append(location)
            elif (gene[key][7][0]<=site<gene[key][8][-1] and (RNAstrand=='.'  or RNAstrand==txstrand)):
                utr3start=gene[key][7]
                utr3end=gene[key][8]
                tx.append(key)
                for j in range(len(utr3start)):
                    if (utr3start[j]<=site<utr3end[j]):
                        if (gene[key][1]==gene[key][2]):
                            location='ncRNA'
                        elif (txstrand=='+'):
                            location='3UTR'
                        elif (txstrand=='-'):
                            location='5UTR'
                        break
                    elif(site>=utr3end[j]):
                        pass
                    elif(site<utr3start[j]):
                        intron=1
                        location='Intron'
                        break
                loctype.append(location)
    
            elif (gene[key][10]<=site<gene[key][11] and (RNAstrand=='.'  or RNAstrand==txstrand)):
                cdsstart=gene[key][1]
                cdsend=gene[key][2]
                tx.append(key)
                for j in range(len(cdsstart)):
                    if (cdsstart[j]<=site<cdsend[j]):
                        location='CDS'
                        length+=site-cdsstart[j]
                        break
                    elif(site>=cdsend[j]):
                        length+=cdsend[j]-cdsstart[j]
                    elif(site<cdsstart[j]):
                        intron=1
                        location='Intron'
                        break
                loctype.append(location)
    locsum=''
    if ('CDS' in loctype):
        locsum='CDS'
    elif ('3UTR' in loctype and '5UTR' in loctype) :
        locsum='5UTR,3UTR'
    elif('5UTR' in loctype):
        locsum='5UTR'
    elif('3UTR' in loctype):
        locsum='3UTR'
    elif('ncRNA' in loctype):
        locsum='ncRNA'
    elif('Intron' in loctype):
        locsum='Intron'
    else:
        pass
    strand_sum=','.join(list(set(strands)))
    
    ref_sense=''
    alt_sense=''
    if (strand_sum=='+' or RNAstrand=='+'):
        ref_sense=ref
        alt_sense=alt
    elif (strand_sum=='-' or RNAstrand=='-'):
        ref_sense=complement[ref]
        alt_sense=complement[alt]
    else:
        pass
    
    siteid=chr+':'+str(site+1)
    
    dbsnp_hit=''
    if (dbsnp.has_key(siteid) and dbsnp[siteid][1]==ref and dbsnp[siteid][2]==alt):
        dbsnp_hit=dbsnp[siteid][0]
    
    radar2_hit='FALSE'
    if (radar2.has_key(siteid) and ref_sense=='A' and alt_sense=='G'):
        radar2_hit='TRUE'
    
    rephit=['']*2
    sitebin2=site/100000
    if (repeat.has_key(chr) and repeat[chr].has_key(sitebin2)):
        for repstart, repend, repname, repclass in repeat[chr][sitebin2]:
            if (repstart<=site<repend ):
                rephit=[repname, repfam]
                break
    else:
        pass
    
    g_count1=a[1].split(',')
    a_count1=a[2].split(',')
    g_count2=a[3].split(',')
    a_count2=a[4].split(',')
    level1=[]
    level2=[]
    total1=0
    total2=0
    num1=0
    num2=0
    for i in range(len(g_count1)):
            g_count1[i]=int(g_count1[i])
            a_count1[i]=int(a_count1[i])
            if (g_count1[i]+a_count1[i]>0):
                level1.append(1.0*g_count1[i]/(g_count1[i]+a_count1[i]))
            else:
                level1.append('NA')  
    for i in range(len(g_count2)):
            g_count2[i]=int(g_count2[i])
            a_count2[i]=int(a_count2[i])
            if (g_count2[i]+a_count2[i]>0):
                level2.append(1.0*g_count2[i]/(g_count2[i]+a_count2[i]))
            else:
                level2.append('NA')
    
    for i in range(len(level1)):
            if (level1[i]!='NA'):
                        total1+=level1[i]
                        num1+=1
                        level1[i]=str(round(level1[i],3))
    for i in range(len(level2)):
            if (level2[i]!='NA'):
                        total2+=level2[i]
                        num2+=1
                        level2[i]=str(round(level2[i],3))
    if (num1>0 and num2>0):
            diff=round(1.0*total1/num1-1.0*total2/num2,3)
    else:
            diff='NA'
            
    onetype=ref_sense+'-'+alt_sense
    if (all_type.has_key(onetype)):
            all_type[onetype]+=1
            fdr=float(a[8])
            #fdr=float(a[8])
            if (fdr<0.05):
                dvr_type[onetype]+=1
                if (radar2_hit=='TRUE'):
                        editing_type[onetype]+=1
                elif (dbsnp_hit!=''):
                        snp_type[onetype]+=1
                else:
                        novel_type[onetype]+=1
    
    if (strand_sum==''):
        strand_sum=RNAstrand
    output.write('\t'.join([chr, str(site+1), ref, alt, RNAstrand, qual]+a[1:5]+a[7:]+[','.join(level1), ','.join(level2), str(diff), ','.join(genename), strand_sum, ref_sense, alt_sense, locsum, dbsnp_hit, radar2_hit, rephit[0], rephit[1]])+'\n')
    
output2=open(args.summary, 'w')
output2.write('\t'.join(['Type (Ref-Alt) on sense strand', 'All variants', 'All DVRs (FDR<0.05)', 'SNP DVRs', 'RNA editing DVRs', 'Novel DVRs'])+'\n')
for onetype in types:
    output2.write('\t'.join([onetype, str(all_type[onetype]), str(dvr_type[onetype]), str(snp_type[onetype]), str(editing_type[onetype]), str(novel_type[onetype])])+'\n')


