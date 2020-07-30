#!/usr/bin/python
import os,sys,time,re
import pandas as pd
aa_123dict = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
              'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}
app='/public/home/sry/DMBS/src/bin/DMBS.sh'
queue='/public/home/sry/DMBS/src/bin/getQ.py'
outdir = '/public/home/sry/DMBS/output'
CAGI_dir = '/public/home/sry/DMBS/dataset/CAGI'
seq_lst_pth = '{CAGI_dir}/seq.lst'.format(CAGI_dir=CAGI_dir)

with open(seq_lst_pth) as f:
    seq_name_lst = [x.strip() for x in f.readlines()]

for seq_name in seq_name_lst:
    fasta_pth = '{CAGI_dir}/fasta/{seq_name}.fasta'.format(CAGI_dir=CAGI_dir,seq_name=seq_name)
    csv_pth = '{CAGI_dir}/{seq_name}.csv'.format(CAGI_dir=CAGI_dir,seq_name=seq_name)
    df = pd.read_csv(csv_pth,low_memory=False)
    print('\n{seq_name}: {_len} mutations'.format(seq_name=seq_name,_len=len(df)))
    sys.stdout.flush()
    for pos in set(df.position.values):
        wt = df.loc[df.position==pos,'start'].values[0]
        mt = [x for x in df.loc[df.position==pos,'end'].values if x in aa_123dict.keys()]
        mut_outdir = '{outdir}/{seq_name}.{pos}{wt}'.format(outdir=outdir,seq_name=seq_name,pos=pos,wt=wt)
        os.makedirs(mut_outdir,exist_ok=True)
        os.system('cp {fasta_pth} {mut_outdir}/seq.fasta'.format(fasta_pth=fasta_pth,mut_outdir=mut_outdir))
        with open('{mut_outdir}/mut_pos.txt'.format(mut_outdir=mut_outdir),'w') as f:
            f.write(str(pos))
        with open('{mut_outdir}/mut_res.txt'.format(mut_outdir=mut_outdir),'w') as f:
            f.write(''.join(set(mt)))

        #qsub
        tag='{seq_name}.{pos}{wt}'.format(seq_name=seq_name,pos=pos,wt=wt)
        walltime = 'walltime=48:00:00'
        qsub_dir = '{mut_outdir}/qsub_master_log'.format(mut_outdir=mut_outdir)
        os.makedirs(qsub_dir,exist_ok=True)
        err='{qsub_dir}/err'.format(qsub_dir=qsub_dir)
        out='{qsub_dir}/out'.format(qsub_dir=qsub_dir)
        run_prog='{qsub_dir}/run_prog.sh'.format(qsub_dir=qsub_dir)

        ## check if the job exists
        user = os.popen("whoami").readline().strip()
        string = '\n'.join(os.popen('/usr/local/bin/qstat -u {user}'.format(user=user)).readlines())
        pattern = re.compile(" " + tag + " ")
        result = pattern.findall(string)
        if result:
            print("job {tag} already running".format(tag=tag))
            sys.stdout.flush()
            continue

        g = open(run_prog, 'w')
        g.writelines('#!/usr/bin/env bash\n')
        g.writelines("echo 'user:' `whoami`\necho 'hostname:' `hostname`\necho 'begin at:' `date`\n")
        g.writelines("echo 'pwd:' `pwd`\n")
        g.writelines('{app} {seq_name}.{pos}{wt} > {qsub_dir}/run.log\n'.format(app=app,seq_name=seq_name,pos=pos,wt=wt,qsub_dir=qsub_dir))
        g.writelines("echo 'end at:' `date`\n")
        g.close()
        
        os.system('chmod +x {run_prog}'.format(run_prog=run_prog))
        os.system(queue)
        os.system('qsub -e {err} -o {out} -l {walltime} -N {tag} {prog}'.format(err=err,out=out,walltime=walltime,tag=tag,prog=run_prog))
        strf = time.strftime("%Y-%m-%d %H:%M:%S")
        print('[{strf}] {run_prog} successfully submitted!'.format(strf=strf, run_prog=run_prog))
        sys.stdout.flush()
        time.sleep(5)