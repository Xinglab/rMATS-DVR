#This script generates the WinBug input files by Simulated RNA-Seq data

import re,os,sys,warnings,numpy,scipy,math,itertools;

from scipy import stats;
from numpy import *;
from multiprocessing import Pool;
from scipy.optimize import fmin_cobyla
from scipy.optimize import fmin_l_bfgs_b
from math import log;

numpy.random.seed(1231);
warnings.filterwarnings('ignore');

#ReadLength
#Dummy length here. Adapt to the new rMATS structure
read_length=50;

#JunctionLength
#Dummy length here. Adapt to the new rMATS structure
junction_length=84;

#Output folder
if len(sys.argv)<3:
       print("Error: Less than two input parameters.");
else:
	output_folder=sys.argv[2];

#MultiProcessor
MultiProcessor=1;
if len(sys.argv)>=4:
	MultiProcessor=int(sys.argv[3]);

#splicing difference cutoff
cutoff=0.1;
if len(sys.argv)>=5:
	cutoff=float(sys.argv[4]);

rho=0.9

#binomial MLE optimization functions
def logit(x):
	if x<0.01:
		x=0.01;
	if x>0.99:
		x=0.99;
	return(log(x/(1-x)));

def logit_list(x_list):
	res=[];
	for x in x_list:
		if x<0.01:
			x=0.01;
		if x>0.99:
			x=0.99;
		res.append(log(x/(1-x)));
	return(res);

#Not use multivar in the MATS LRT 
def myfunc_multivar(x,*args):
	psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
	sum1=0;sum2=0;
	for i in range(len(psi1)):
		sum1+=pow(logit((psi1[i]-psi2[i])/2+0.5)-logit((x[0]-x[1])/2+0.5),2)
	sum1=sum1/var2/2;
	return(sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x[0]),2)+pow(stats.norm.ppf(x[1]),2)-2*rho*stats.norm.ppf(x[0])*stats.norm.ppf(x[1])));
	#return(sum1);

#Not use multivar in the MATS LRT 
def myfunc_multivar_der(x,*args):
	psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
	sum1=0;sum2=0;
	for i in range(len(psi1)):
		sum1+=-2*(logit((psi1[i]-psi2[i])/2+0.5)-logit((x[0]-x[1])/2+0.5))/((x[0]-x[1])/2+0.5)/(0.5-(x[0]-x[1])/2)*0.5;
	sum1=sum1/var1/2;
	res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[0])-2*rho*stats.norm.ppf(x[1]))/stats.norm.pdf(stats.norm.ppf(x[0]));
	for i in range(len(psi2)):
		sum2+=-2*(logit((psi1[i]-psi2[i])/2+0.5)-logit((x[0]-x[1])/2+0.5))/((x[0]-x[1])/2+0.5)/(0.5-(x[0]-x[1])/2)*(-0.5);
	sum2=sum2/var2/2;
	res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[1])-2*rho*stats.norm.ppf(x[0]))/stats.norm.pdf(stats.norm.ppf(x[1]));
	return(numpy.array([res1,res2]));

#note to me: change this in unpaired rMATS
def myfunc_1(x,*args):
	I1=args[0][0];I2=args[0][1];S1=args[1][0];S2=args[1][1];beta1=args[2][0];beta2=args[2][1];var=args[3];effective_inclusion_length=args[4];effective_skipping_length=args[5];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	new_psi1=inclusion_length*(x+cutoff)/(inclusion_length*(x+cutoff)+skipping_length*(1-(x+cutoff)));new_psi2=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x));
	binomial_sum=-1*(I1*log(new_psi1)+S1*log(1-new_psi1)+I2*log(new_psi2)+S2*log(1-new_psi2));
	multivar_sum=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x+cutoff),2)+pow(stats.norm.ppf(x),2)-2*rho*stats.norm.ppf(x+cutoff)*stats.norm.ppf(x))
	return(binomial_sum+multivar_sum);

def myfunc_der_1(x,*args):
	I1=args[0][0];I2=args[0][1];S1=args[1][0];S2=args[1][1];beta1=args[2][0];beta2=args[2][1];var=args[3];effective_inclusion_length=args[4];effective_skipping_length=args[5];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	new_psi1=inclusion_length*(x+cutoff)/(inclusion_length*(x+cutoff)+skipping_length*(1-(x+cutoff)));new_psi2=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x));
	new_psi1_der=inclusion_length*skipping_length/pow(inclusion_length*(x+cutoff)+skipping_length*(1-(x+cutoff)),2);
	new_psi2_der=inclusion_length*skipping_length/pow(inclusion_length*x+skipping_length*(1-x),2);
	res1=-1*(I1/new_psi1-S1/(1-new_psi1))*new_psi1_der;
	res1+=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x+cutoff)-2*rho*stats.norm.ppf(x))/stats.norm.pdf(stats.norm.ppf(x+cutoff))
	res2=-1*(I2/new_psi2-S2/(1-new_psi2))*new_psi2_der;
	res2+=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x)-2*rho*stats.norm.ppf(x+cutoff))/stats.norm.pdf(stats.norm.ppf(x));
	return(numpy.array(res1+res2));

def myfunc_2(x, *args):
	I1=args[0][0];I2=args[0][1];S1=args[1][0];S2=args[1][1];beta1=args[2][0];beta2=args[2][1];var=args[3];effective_inclusion_length=args[4];effective_skipping_length=args[5];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	new_psi1=inclusion_length*(x)/(inclusion_length*(x)+skipping_length*(1-(x)));new_psi2=inclusion_length*(x+cutoff)/(inclusion_length*(x+cutoff)+skipping_length*(1-(x+cutoff)));
	binomial_sum=-1*(I1*log(new_psi1)+S1*log(1-new_psi1)+I2*log(new_psi2)+S2*log(1-new_psi2));
	multivar_sum=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x),2)+pow(stats.norm.ppf(x+cutoff),2)-2*rho*stats.norm.ppf(x)*stats.norm.ppf(x+cutoff))
	return(binomial_sum+multivar_sum);

def myfunc_der_2(x,*args):
	I1=args[0][0];I2=args[0][1];S1=args[1][0];S2=args[1][1];beta1=args[2][0];beta2=args[2][1];var=args[3];effective_inclusion_length=args[4];effective_skipping_length=args[5];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	new_psi1=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x));new_psi2=inclusion_length*(x+cutoff)/(inclusion_length*(x+cutoff)+skipping_length*(1-(x+cutoff)));
	new_psi1_der=inclusion_length*skipping_length/pow(inclusion_length*x+skipping_length*(1-x),2);
	new_psi2_der=inclusion_length*skipping_length/pow(inclusion_length*(x+cutoff)+skipping_length*(1-(x+cutoff)),2);
	res1=-1*(I1/new_psi1-S1/(1-new_psi1))*new_psi1_der;
	res1+=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x)-2*rho*stats.norm.ppf((x+cutoff)))/stats.norm.pdf(stats.norm.ppf(x))
	res2=-1*(I2/new_psi2-S2/(1-new_psi2))*new_psi2_der;
	res2+=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf((x+cutoff))-2*rho*stats.norm.ppf(x))/stats.norm.pdf(stats.norm.ppf((x+cutoff)));
	return(numpy.array(res1+res2));

def myfunc_individual(x,*args):
	I1=args[0][0];I2=args[0][1];S1=args[1][0];S2=args[1][1];beta1=args[2][0];beta2=args[2][1];var=args[3];effective_inclusion_length=args[4];effective_skipping_length=args[5];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	new_psi1=inclusion_length*x[0]/(inclusion_length*x[0]+skipping_length*(1-x[0]));new_psi2=inclusion_length*x[1]/(inclusion_length*x[1]+skipping_length*(1-x[1]));
	binomial_sum=-1*(I1*log(new_psi1)+S1*log(1-new_psi1)+I2*log(new_psi2)+S2*log(1-new_psi2));
	multivar_sum=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x[0]),2)+pow(stats.norm.ppf(x[1]),2)-2*rho*stats.norm.ppf(x[0])*stats.norm.ppf(x[1]))
	return(binomial_sum+multivar_sum);

def myfunc_individual_der(x,*args):
	I1=args[0][0];I2=args[0][1];S1=args[1][0];S2=args[1][1];beta1=args[2][0];beta2=args[2][1];var=args[3];effective_inclusion_length=args[4];effective_skipping_length=args[5];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	new_psi1=inclusion_length*x[0]/(inclusion_length*x[0]+skipping_length*(1-x[0]));new_psi2=inclusion_length*x[1]/(inclusion_length*x[1]+skipping_length*(1-x[1]));
	new_psi1_der=inclusion_length*skipping_length/pow(inclusion_length*x[0]+skipping_length*(1-x[0]),2);
	new_psi2_der=inclusion_length*skipping_length/pow(inclusion_length*x[1]+skipping_length*(1-x[1]),2);
	res1=-1*(I1/new_psi1-S1/(1-new_psi1))*new_psi1_der;
	res1+=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[0])-2*rho*stats.norm.ppf(x[1]))/stats.norm.pdf(stats.norm.ppf(x[0]))
	res2=-1*(I2/new_psi2-S2/(1-new_psi2))*new_psi2_der;
	res2+=0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[1])-2*rho*stats.norm.ppf(x[0]))/stats.norm.pdf(stats.norm.ppf(x[1]));
	return(numpy.array([res1,res2]));

def myfunc_likelihood(x, args):
	I1=args[0][0];I2=args[0][1];S1=args[1][0];S2=args[1][1];beta1=args[2][0];beta2=args[2][1];var=args[3];
	sum=0;N1=I1+S1;N2=I2+S2;
	if (N1+N2)==0:
		return(0);
	sum+=-0.5*((I1-N1*x[0])*(I1-N1*x[0])/(N1*x[0])+(S1-N1*(1-x[0]))*(S1-N1*(1-x[0]))/(N1*(1-x[0])));
	sum+=-0.5*((I2-N2*x[1])*(I2-N2*x[1])/(N2*x[1])+(S2-N2*(1-x[1]))*(S2-N2*(1-x[1]))/(N2*(1-x[1])));
	sum+=pow(logit(beta1)-logit(beta2)-logit(x[0])+logit(x[1]),2);
	return(sum);

def MLE_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
	psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
	while((iter_cutoff>0.01)&(count<=iter_maxrun)):
		count+=1;
		#iteration of beta
		beta_0=sum(psi1)/len(psi1);
		beta_1=sum(psi2)/len(psi2);
		var1=0;var2=0;
		current_sum=0;likelihood_sum=0;
		new_psi1=[];new_psi2=[];
		if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
			xopt = fmin_l_bfgs_b(myfunc_1,[sum(psi2)/len(psi2)],myfunc_der_1,args=[[i1[0],i2[0]],[s1[0],s2[0]],[beta_0,beta_1],var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.001,0.999-cutoff]],iprint=-1)
			theta2 = max(min(float(xopt[0]),1-cutoff),0);theta1=theta2+cutoff;
		else:#minize psi1 if this is the case
			xopt = fmin_l_bfgs_b(myfunc_2,[sum(psi1)/len(psi1)],myfunc_der_2,args=[[i1[0],i2[0]],[s1[0],s2[0]],[beta_0,beta_1],var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.001,0.999-cutoff]],iprint=-1)
			theta1 =	 max(min(float(xopt[0]),1-cutoff),0);theta2=theta1+cutoff;
		#Debug;print('constrain_1xopt');print('theta');print(theta1);print(theta2);print(xopt);
		current_sum+=float(xopt[1]);
		new_psi1.append(theta1);new_psi2.append(theta2);
		psi1=new_psi1;psi2=new_psi2;
		if count>1:
			iter_cutoff=abs(previous_sum-current_sum)/abs(previous_sum);
		previous_sum=current_sum;
	#Debug;print('constrain');print(theta1);print(theta2);print(psi1);print(psi2);print(current_sum);print(likelihood_sum);
	print('constrain');print(xopt);print(theta1);print(theta2);
	return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);

def MLE_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
	psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
	while((iter_cutoff>0.01)&(count<=iter_maxrun)):
		count+=1;
		#iteration of beta
		beta_0=sum(psi1)/len(psi1);
		beta_1=sum(psi2)/len(psi2);
		var1=0;var2=0;
		current_sum=0;likelihood_sum=0;
		new_psi1=[];new_psi2=[];
		#Debug;print('unconstrain_1xopt');
		for i in range(len(psi1)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i],psi2[i]],myfunc_individual_der,args=[[i1[i],i2[i]],[s1[i],s2[i]],[beta_0,beta_1],var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99],[0.01,0.99]],iprint=-1);
			new_psi1.append(float(xopt[0][0]));current_sum+=float(xopt[1]);
			new_psi2.append(float(xopt[0][1]));
			#Debug;print(xopt);
			likelihood_sum+=myfunc_likelihood([new_psi1[i],new_psi2[i]],[[i1[i],i2[i]],[s1[i],s2[i]],[beta_0,beta_1],var1]);
		psi1=new_psi1;psi2=new_psi2;
		#Debug;print('count');print(count);print('previous_sum');print(previous_sum);print('current_sum');print(current_sum);
		if count>1:
			iter_cutoff=abs(previous_sum-current_sum)/abs(previous_sum);
		previous_sum=current_sum;
	if count>iter_maxrun:
		return([current_sum,[psi1,psi2,0,0,var1,var2]]);
	print('unconstrain');print(xopt);
	return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);

#Random Sampling Function
def likelihood_test(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length,flag):
	if flag==0:
		return(1);
	else:
		res=MLE_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length);
		if abs(res[1][2]-res[1][3])<=cutoff:
			#Debug;print('1<=cutoff');print(res);print((res[1][2]-res[1][3]));
			return(1);
		else:
			res_constrain=MLE_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length);
			#Debug;print('2>cutoff');print('res');print(res);print('res_constrain');print(res_constrain);
			#Debug;print(abs(res_constrain[0]-res[0]));print('2end');
			return(1-scipy.stats.chi2.cdf(2*(abs(res_constrain[0]-res[0])),1));

#MultiProcessorFunction
def MultiProcessorPool(n_original_diff):
	i1=n_original_diff[0];i2=n_original_diff[1];s1=n_original_diff[2];s2=n_original_diff[3];effective_inclusion_length=n_original_diff[4];effective_skipping_length=n_original_diff[5];flag=n_original_diff[6];
	P=likelihood_test(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length,flag);
	return(P);

#Function for vector handling
def vec2float(vec):
	res=[];
	for i in vec:
		res.append(float(i));
	return(res);

def vec2sum(vec):
	res=0;
	for i in vec:
		res+=float(i);
	return([res]);
	
def vecprod(vec):
	res=1;
	for i in vec:
		res=res*i;
	return(res);

def vec2remove0psi(inc,skp):
	res1=[];res2=[];
	for i in range(len(inc)):
		if (inc[i]!=0) | (skp[i]!=0):
			res1.append(inc[i]);res2.append(skp[i]);
	return([res1,res2]);

def vec2psi(inc,skp,effective_inclusion_length,effective_skipping_length):
	psi=[];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	for i in range(len(inc)):
		if (float(inc[i])+float(skp[i]))==0:
			psi.append(0.5);
		else:
			psi.append(float(inc[i])/inclusion_length/(float(inc[i])/inclusion_length+float(skp[i])/skipping_length));
	return(psi);

def vec210(vec):
	res=[];
	for i in vec:
		if i>0:
			res.append(1);
		else:
			res.append(-1);
	return(res);

def myttest(vec1,vec2):
	if (len(vec1)==1) & (len(vec2)==1):
		res=stats.ttest_ind([vec1[0],vec1[0]],[vec2[0],vec2[0]]);
	else:
		res=stats.ttest_ind(vec1,vec2);
	return(res);

ifile=open(sys.argv[1]);
title=ifile.readline();
#analyze the title of the inputed data file to find the information of how much simulation are involved
#the min simulated round is 10, each time it increases by 10 times
element=re.findall('[^ \t\n]+',title);
ofile=open(output_folder+'/rMATS_Result_P.txt','w');
ofile.write(title[:-1]+'\tPValue'+'\n');

list_n_original_diff=[];probability=[];psi_list_1=[];psi_list_2=[];rho_list=[];psi1_for_rho_list=[];psi2_for_rho_list=[];
ilines=ifile.readlines();
for i in range(len(ilines)):
	element=re.findall('[^ \t\n]+',ilines[i]);
	inc1=re.findall('[^,]+',element[1]);skp1=re.findall('[^,]+',element[2]);inc2=re.findall('[^,]+',element[3]);skp2=re.findall('[^,]+',element[4]);
	#Dummy effective_inclusion_length and flanking exon length here. Adapt to the new rMATS structure
	effective_inclusion_length=int(element[5]);
	effective_skipping_length=int(element[6]);
	#inc1=vec2float(inc1);skp1=vec2float(skp1);inc2=vec2float(inc2);skp2=vec2float(skp2);
	inc1=vec2sum(inc1);skp1=vec2sum(skp1);inc2=vec2sum(inc2);skp2=vec2sum(skp2);
	if ((vecprod(inc1)+vecprod(skp1))==0) | ((vecprod(inc2)+vecprod(skp2))==0):
		list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0]);
	else:
		psi1=vec2psi(inc1,skp1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(inc2,skp2,effective_inclusion_length,effective_skipping_length);
		for i in range(len(psi1)):
			if len(psi1_for_rho_list)<=i: 
				psi1_for_rho_list.append([]);
			psi1_for_rho_list[i].append(psi1[i]);
		for i in range(len(psi2)):
			if len(psi2_for_rho_list)<=i:
				psi2_for_rho_list.append([]);
			psi2_for_rho_list[i].append(psi2[i]);
		psi_list_1.append(sum(inc1)/(sum(inc1)+sum(skp1)));
		psi_list_2.append(sum(inc2)/(sum(inc2)+sum(skp2)));
		#temp1=vec2remove0psi(inc1,skp1);temp2=vec2remove0psi(inc2,skp2);
		#inc1=temp1[0];skp1=temp1[1];inc2=temp2[0];skp2=temp2[1];
		list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,1]);
	#if i>2:
	#	break;

#rho_list for paired data
for i in range(len(psi1_for_rho_list)):
	this_rho=stats.pearsonr(numpy.array(psi1_for_rho_list[i]),numpy.array(psi2_for_rho_list[i]));this_rho=this_rho[0];
	if this_rho>0.9:
		this_rho=0.9;
	rho_list.append(this_rho);

rho=stats.pearsonr(numpy.array(psi_list_1),numpy.array(psi_list_2));rho=rho[0];
if rho>0.9:
	rho=0.9;

#rho_list=[0.95,0.95,0.95,0.95];
#rho_list=[0,0,0,0];
rho=0.9;
#print('rho');print(rho);

if MultiProcessor>1:
	pool=Pool(processes=MultiProcessor);
	probability=pool.map(MultiProcessorPool,list_n_original_diff);
else:
	for i in range(len(list_n_original_diff)):
		#print(list_n_original_diff[i]);
		probability.append(MultiProcessorPool(list_n_original_diff[i]));

#print(probability);
index=0;
for i in range(len(ilines)):
    element=re.findall('[^ \t\n]+',ilines[i]);
    ofile.write(ilines[i][:-1]+'\t'+str(probability[i])+'\n');
ofile.close();
