#!/usr/bin/python
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, type=str,help='inputfile')
parser.add_argument('-P', '--outpssm', required=True, type=str,help='output_pssmfile')
parser.add_argument('-H', '--outhhw', required=True, type=str,help='output_hhwfile')
args=parser.parse_args()

inputfile       = args.input
pssm_outputfile = args.outpssm
hhw_outputfile  = args.outhhw


PSSM_col=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
f=open(inputfile,"r")
a=f.readlines()
query=a[0].strip()
pd_query=pd.DataFrame(list(query),columns=["Res"])
Num=len(a)
L=len(query)
num_count=np.zeros((L,20))
for lines in a:
	lines=lines.strip()
	for i,r in enumerate(lines):
		if r in PSSM_col:
			j=PSSM_col.index(r)
			num_count[i][j]=num_count[i][j]+1

f.close()

N=[]
N_all=[]
for x in num_count:
	n=np.sum(x!=0)
	n_all=np.sum(x)
	N.append(n)
	N_all.append(int(n_all))

PSSM_num=np.dot(np.diag(100.0/np.array(N_all)),num_count).astype("int")

pd_PSSM=pd.DataFrame(PSSM_num,columns=PSSM_col)
pd_PSSM=pd.concat([pd_query,pd_PSSM],axis=1)
pd_PSSM.to_csv(pssm_outputfile)

hhw=np.zeros((Num,L))
for i,lines in enumerate(a):
	lines=lines.strip()
	for j,r in enumerate(lines):
		if r in PSSM_col:
			k=PSSM_col.index(r)
			w_jk=num_count[j,k]
			w_j=N[j]
			hhw[i,j]=1.0/(w_jk*w_j)

sum_N=np.sum(hhw,axis=1)
wt_N=sum_N/np.sum(sum_N)

PSSM_hhw=np.zeros((L,20))
for i,lines in enumerate(a):
	lines=lines.strip()
	for j,r in enumerate(lines):
		if r in PSSM_col:
			k=PSSM_col.index(r)
			PSSM_hhw[j][k]=PSSM_hhw[j][k]+wt_N[i]

wt_all=np.sum(PSSM_hhw,axis=1)
PSSM_hhw=np.dot(np.diag(100.0/np.array(wt_all)),PSSM_hhw).astype("int")
pd_hhw=pd.DataFrame(PSSM_hhw,columns=PSSM_col)
pd_hhw=pd.concat([pd_query,pd_hhw],axis=1)
pd_hhw.to_csv(hhw_outputfile)
