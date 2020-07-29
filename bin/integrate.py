#!/usr/bin/python
import os, sys, csv, re
import numpy as np
import pandas as pd
from scipy import stats

mut_pos=int(sys.argv[1].strip())
mut_res=sys.argv[2].strip()
indir=sys.argv[3].strip()
blosum62=sys.argv[4].strip()

outdir="{indir}/{mut_res}".format(indir=indir, mut_res=mut_res)
logdir="{outdir}/feature_log".format(outdir=outdir)
os.makedirs(logdir, exist_ok=True)

B=pd.read_csv(blosum62,header=0,sep=',',index_col=0)

#tag="PSSM_"+name
fasta_file_name=indir+"/seq.fasta"

#read the psi-blast out_ascii_pssm file into a pssm list
#   the list: pos stores the position information (range(n) where n is the length of the query)
#   the list: S stores the residues in query
#   the list: log_odds stores the query's MSA pssm matrix which has changed into log_odds type
#   the list: frequence stores the query's MSA pssm matrix which is in frequence type
def read_pssm(filename):
    f=pd.read_csv(filename)
    l=len(f)
    keys=list(f.columns[2:])
    pos=np.array(range(l))+1
    S=list(f["Res"])
    frequence=np.array(f[keys],dtype="int")
    return [pos,S,keys,frequence]

def read_SVM(filename):
    f=open(filename,"r")
    a=f.readlines()
    f.close()
    b=a[0].strip().split()[1]
    target_column=1 if b=='1' else 2
    Bp=[]
    qua=[]
    for i,x in enumerate(a[1:]):
        if x[0]=='1':
            q=x.strip().split()[target_column]
            Bp.append(i+1)
            qua.append(q)
    if Bp==[]:
        qua=max(map(lambda x:x.strip().split()[target_column],a[1:]))
    return [Bp,qua]

def read_ssite(filename):
    f=open(filename)
    a=f.readlines()
    f.close()
    l=len(a)
    qua=[]
    Bp=[]
    for x in a:
        slice=x.strip().split('\t')
        q=float(slice[0])
        pos=slice[-2].split(',')
        for d in range(len(pos)-1):
            ps=int(pos[d])
            pe=int(pos[d+1])
            length=pe-ps
            if length <= 5:
                Bp.extend(range(ps,pe))
                qua.extend([q]*length)
            else:
                Bp.append(ps)
                qua.append(q)
        qua.append(q)
        Bp.append(int(pos[-1]))
    return [Bp,qua]


#calculate the BLOSUM62 score
def BLO_score(B,keys,S,pssm_row,hhw=False):
    score=0
    for i,res in enumerate(keys):
        score+= B[res][S]*pssm_row[i]
    return score/100.0

#calculate the entropy score
def data_entropy(pssm_row):
    pssm_row=np.array(pssm_row,dtype="float")
    p_data=pssm_row/sum(pssm_row)
    ent=stats.entropy(p_data)
    return ent

def KL(P,Q):
    """ Epsilon is used here to avoid conditional code for
    checking that neither P nor Q is equal to 0. """
    epsilon = 0.00001
    
    #You may want to instead make copies to avoid changing the np arrays.
    P = P+epsilon
    Q = Q+epsilon
    
    divergence = np.sum(P*np.log(P/Q))
    return divergence


def cal_index(Bp,qua,pssm_ini,pssm_mut,keys_ini=None,sequence=None,mut_pos=None,mut_res=None,kind='kl'):
    p_max=0
    tar_max=0
    q=min(qua)
    for k,pos in enumerate(Bp):
        pssm_row_ini=np.array(pssm_ini[pos-1],dtype="float")
        pssm_row_mut=np.array(pssm_mut[pos-1],dtype="float")
        if kind=='kl':
            tar_=KL(pssm_row_ini,pssm_row_mut)
        elif kind=='entropy':
            tar_=abs(data_entropy(pssm_row_ini)-data_entropy(pssm_row_mut))
        elif kind=="BLOSUM":
            if pos == mut_pos:
                tar_=BLO_score(B,keys_ini,sequence[pos-1],pssm_ini[pos-1],hhw=False)-BLO_score(B,keys_ini,mut_res,pssm_mut[pos-1],hhw=False)
            else:
                tar_=BLO_score(B,keys_ini,sequence[pos-1],pssm_ini[pos-1],hhw=False)-BLO_score(B,keys_ini,sequence[pos-1],pssm_mut[pos-1],hhw=False)
            tar_=abs(tar_)
        else:
            print("Not correct objective function!")
        if tar_ > tar_max:
            tar_max=tar_
            p_max=pos-1
            q=qua[k]
    return [tar_max,p_max,q]



def data_slice(Ent,FQ,keys,pos,ini,mut):
    ent_ini=Ent[0]
    ent_mut=Ent[1]
    ent_hhw_ini=Ent[2]
    ent_hhw_mut=Ent[3]
    fq_ini=FQ[0]
    fq_mut=FQ[1]
    ini_hhw=FQ[2]
    mut_hhw=FQ[3]
    keys_ini=keys[0]
    keys_hhw=keys[1]
    entropy_sw_ini=ent_ini[pos]
    entropy_sw_mut=ent_mut[pos]
    entropy_sw_delta=entropy_sw_ini-entropy_sw_mut
    entropy_hhw_ini=ent_hhw_ini[pos]
    entropy_hhw_mut=ent_hhw_mut[pos]
    entropy_hhw_delta=entropy_hhw_ini-entropy_hhw_mut
    
    BLOSUM62_sw_ini=BLO_score(B,keys_ini,ini,fq_ini[pos],hhw=False)
    BLOSUM62_sw_mut=BLO_score(B,keys_ini,mut,fq_mut[pos],hhw=False)
    BLOSUM62_sw_delta=BLOSUM62_sw_ini-BLOSUM62_sw_mut
    BLOSUM62_hhw_ini=BLO_score(B,keys_hhw,ini,ini_hhw[pos],hhw=True)
    BLOSUM62_hhw_mut=BLO_score(B,keys_hhw,mut,mut_hhw[pos],hhw=True)
    BLOSUM62_hhw_delta=BLOSUM62_hhw_ini-BLOSUM62_hhw_mut
    
    mean_ini=np.mean(np.array(ent_ini))
    std_ini=np.std(np.array(ent_ini))
    Z_score_ini=(entropy_sw_ini-mean_ini)/std_ini
    
    mean_mut=np.mean(np.array(ent_mut))
    std_mut=np.std(np.array(ent_mut),ddof=1)

    Z_score_mut=(entropy_sw_mut-mean_mut)/std_mut
    
    mean_delta=mean_ini-mean_mut
    std_delta=std_ini-std_mut
    Z_score_delta=Z_score_ini-Z_score_mut
    
    d_slice=[entropy_sw_ini,entropy_sw_mut,entropy_sw_delta,
        entropy_hhw_ini,entropy_hhw_mut,entropy_sw_delta,
        BLOSUM62_sw_ini,BLOSUM62_sw_mut,BLOSUM62_sw_delta,
        BLOSUM62_hhw_ini,BLOSUM62_hhw_mut,BLOSUM62_hhw_delta,
        mean_ini,std_ini,Z_score_ini,
        mean_mut,std_mut,Z_score_mut,
        mean_delta,std_delta,Z_score_delta]
    return d_slice

a=["bla","hhm"]
data=[]
thr=0
ws=2*thr+1
len_col=21
len_slice=ws*len_col
g=open(fasta_file_name,"r")
b=g.readlines()
g.close()
sequence=""
for lines in b:
    if re.match(r"\w",lines):
        sequence=sequence+lines.strip()

for types in a:
    #We can find the pssm files in the file whose path is path_to_file
    try:
        [pos_ini,S_ini,keys_ini,fq_ini]=read_pssm(indir+'/seq_{tp}.pssm'.format(tp=types))
        [pos_mut,S_mut,keys_mut,fq_mut]=read_pssm(outdir+'/mut_{tp}.pssm'.format(tp=types))
        [_,_,keys_hhw,ini_hhw]=read_pssm(indir+'/seq_{tp}.hhw'.format(tp=types))
        [_,_,keys_hhw,mut_hhw]=read_pssm(outdir+'/mut_{tp}.hhw'.format(tp=types))
        len_seq=len(pos_ini)-1
    except IOError:
        f=open(logdir+"/warnings.info","w+")
        f.write("Error in reading pssm")
        f.close()
        new_data_list=[]
        d_slice=list([0])*len_slice
        new_data_list.extend(d_slice)
    else:
        ent_ini=[]
        ent_mut=[]
        for k in range(len(S_ini)):
            ent_ini.append(data_entropy(fq_ini[k]))
            ent_mut.append(data_entropy(fq_mut[k]))
        
        ent_hhw_ini=[]
        ent_hhw_mut=[]
        for k in range(len(S_ini)):
            ent_hhw_ini.append(data_entropy(ini_hhw[k]))
            ent_hhw_mut.append(data_entropy(mut_hhw[k]))
        
        Ent=list([ent_ini,ent_mut,ent_hhw_ini,ent_hhw_mut])
        FQ=list([fq_ini,fq_mut,ini_hhw,mut_hhw])
        keys=list([keys_ini[:20],keys_hhw])
        new_data_list=[]
        for x in range(ws):
            pos=int(mut_pos)-1-thr+x
            if (pos < 0) or (pos > len_seq):
                d_slice=list([0])*len_col
            else:
                if pos==(int(mut_pos)-1):
                    ini=sequence[pos]
                    mut=mut_res
                else:
                    ini=sequence[pos]
                    mut=ini
                try:
                    d_slice=data_slice(Ent,FQ,keys,pos,ini,mut)
                except KeyError as e:
                    print(e)
                    d_slice=list([0])*len_col
                else:
                    pass
            new_data_list.extend(d_slice)
        
    data.extend(new_data_list)


##################################################    pssm added         
"""
#############################################################
len_col=9
len_slice=len_col
ssite_file_name=outdir+"/Bsites_fpt.dat"
[Bp,qua]=read_ssite(ssite_file_name)
try:
    [Bp,qua]=read_ssite(ssite_file_name)
    len_seq=len(pos_ini)-1
except IOError:
    print ("No file in "+outdir)
    new_data_list=[]
    d_slice=list([0])*len_slice
    new_data_list.extend(d_slice)
else:
    [t_kl,pos_max_kl,q_kl]=cal_index(Bp,qua,fq_ini,fq_mut,kind='kl')
    [t_ent,pos_max_ent,q_kl]=cal_index(Bp,qua,fq_ini,fq_mut,kind='entropy')
    [t_B,pos_max_B,q_B]=cal_index(Bp,qua,fq_ini,fq_mut,keys_ini=keys_ini,sequence=sequence,mut_pos=pos,mut_res=mut_res,kind="BLOSUM")
    new_data_list=[]
    d_slice=[t_kl,pos_max_kl,q_kl,t_ent,pos_max_ent,q_kl,t_B,pos_max_B,q_B]
    new_data_list.extend(d_slice)

data.extend(new_data_list)

#################################################     ssite bind added  ############################################################
if (mut_pos in Bp):
    data.append(0)
else:
    d=int(min(abs(np.array(Bp)-mut_pos)))
    data.append(d)
##################################################    binding site distance added     #########################################################
len_col=18
len_slice=len_col
DNA_ssite_file_name=outdir+"/5D-out"
RNA_ssite_file_name=outdir+"/5R-out"
try:
    [Bp_D,qua_D]=read_SVM(DNA_ssite_file_name)
    [Bp_R,qua_R]=read_SVM(RNA_ssite_file_name)
    len_seq=len(pos_ini)-1
except IOError:
    print ("No file in "+outdir)
    new_data_list=[]
    d_slice=list([0])*len_slice
    new_data_list.extend(d_slice)
else:
    def gs(Bp,qua):
        [t_kl,pos_max_kl,q_kl]=cal_index(Bp,qua,fq_ini,fq_mut,kind='kl')
        [t_ent,pos_max_ent,q_kl]=cal_index(Bp,qua,fq_ini,fq_mut,kind='entropy')
        [t_B,pos_max_B,q_B]=cal_index(Bp,qua,fq_ini,fq_mut,keys_ini=keys_ini,sequence=sequence,mut_pos=pos,mut_res=mut_res,kind="BLOSUM")
        d_slice=[t_kl,pos_max_kl,q_kl,t_ent,pos_max_ent,q_kl,t_B,pos_max_B,q_B]
        return d_slice
    
    d_slice_D=gs(Bp_D,qua_D)
    d_slice_R=gs(Bp_R,qua_R)
    
    new_data_list=[]
    new_data_list.extend(d_slice_D)
    new_data_list.extend(d_slice_R)

data.extend(new_data_list)
####################################################    svmsite bind added   ############################################################
try:
    f=open("%s/5D-out"%outdir,"r")
except IOError:
    print("NO file in %s"%outdir)
    D_con=-1
    D_pro=0
    R_con=-1
    R_pro=0
else:
    a=f.readlines()
    f.close()
    info=a[mut_pos].split()
    D_con=info[0]
    D_pro=info[-1]
    f=open("%s/5R-out"%outdir,"r")
    a=f.readlines()
    f.close()
    info=a[mut_pos].split()
    R_con=info[0]
    R_pro=info[1]

data.extend([D_con,D_pro,R_con,R_pro])
"""

#labels=["etp_eql_ini_bla","etp_eql_mut_bla","etp_eql_delta_bla","etp_hhw_ini_bla","etp_hhw_mut_bla","etp_hhw_delta_bla","blo_eql_ini_bla","blo_eql_mut_bla","blo_eql_delta_bla","blo_hhw_ini_bla","blo_hhw_mut_bla","blo_hhw_delta_bla","mean_ini_bla","std_ini_bla","Z_ini_bla","mean_mut_bla","std_mut_bla","Z_mut_bla","mean_delta_bla","std_delta_bla","Z_delta_bla","etp_eql_ini_0_hhm","etp_eql_mut_0_hhm","etp_eql_delta_0_hhm","etp_hhw_ini_0_hhm","etp_hhw_mut_0_hhm","etp_hhw_delta_0_hhm","blo_eql_ini_0_hhm","blo_eql_mut_0_hhm","blo_eql_delta_0_hhm","blo_hhw_ini_0_hhm","blo_hhw_mut_0_hhm","blo_hhw_delta_0_hhm","mean_ini_0_hhm","std_ini_0_hhm","Z_ini_0_hhm","mean_mut_0_hhm","std_mut_0_hhm","Z_mut_0_hhm","mean_delta_0_hhm","std_delta_0_hhm","Z_delta_0_hhm",
#"kl_max_ssite","kl_pos_ssite","kl_qua_ssite","ent_max_ssite","ent_pos_ssite","ent_qua_ssite","B_max_ssite","B_pos_ssite","B_qua_ssite","Bind",
#"kl_max_d_svm","kl_pos_d_svm","kl_qua_d_svm","ent_max_d_svm","ent_pos_d_svm","ent_qua_d_svm","B_max_d_svm","B_pos_d_svm","B_qua_d_svm","kl_max_r_svm","kl_pos_r_svm","kl_qua_r_svm","ent_max_r_svm","ent_pos_r_svm","ent_qua_r_svm","B_max_r_svm","B_pos_r_svm","B_qua_r_svm",
#"D_con","D_pro","R_con","R_pro"]
labels=["etp_eql_ini_bla","etp_eql_mut_bla","etp_eql_delta_bla","etp_hhw_ini_bla","etp_hhw_mut_bla","etp_hhw_delta_bla","blo_eql_ini_bla","blo_eql_mut_bla","blo_eql_delta_bla","blo_hhw_ini_bla","blo_hhw_mut_bla","blo_hhw_delta_bla","mean_ini_bla","std_ini_bla","Z_ini_bla","mean_mut_bla","std_mut_bla","Z_mut_bla","mean_delta_bla","std_delta_bla","Z_delta_bla","etp_eql_ini_0_hhm","etp_eql_mut_0_hhm","etp_eql_delta_0_hhm","etp_hhw_ini_0_hhm","etp_hhw_mut_0_hhm","etp_hhw_delta_0_hhm","blo_eql_ini_0_hhm","blo_eql_mut_0_hhm","blo_eql_delta_0_hhm","blo_hhw_ini_0_hhm","blo_hhw_mut_0_hhm","blo_hhw_delta_0_hhm","mean_ini_0_hhm","std_ini_0_hhm","Z_ini_0_hhm","mean_mut_0_hhm","std_mut_0_hhm","Z_mut_0_hhm","mean_delta_0_hhm","std_delta_0_hhm","Z_delta_0_hhm"]
write_data=[labels,data]

with open(outdir+'/DMBS_feature.csv','w+') as csvfile:
    writer=csv.writer(csvfile)
    writer.writerows(write_data)
