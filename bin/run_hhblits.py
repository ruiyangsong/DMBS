#!/usr/bin/python

import os, sys, time, re

hhblits_app, libdir, seq_path, outdir, data_name, out_a3m, out_hhm, queue,\
out_aln,\
convert_aln_app, tmp_a3m,\
hhfilter_app, filtered_a3m, filtered_aln,\
pssm_app, pssm_a3m, hhw_a3m,\
mod_pth = sys.argv[1:] = sys.argv[1:]

logdir = '{outdir}/qsub_hhblits_log'.format(outdir=outdir)
os.makedirs(logdir, exist_ok=True)

tag      = "hh.{name}".format(name=data_name)
err      = "{logdir}/err".format(logdir=logdir)
out      = "{logdir}/out".format(logdir=logdir)
walltime = "walltime=05:00:00"
run_prog = "{logdir}/run_prog.sh".format(logdir=logdir)


## check if the job exists
user = os.popen("whoami").readline().strip()
string = '\n'.join(os.popen('/usr/local/bin/qstat -f | grep -B 1 "Job_Owner.*{user}"'.format(user=user)).readlines())
pattern = re.compile(tag)
result = pattern.findall(string)
if result:
    print("job {tag} already running".format(tag=tag))
    exit(0)


f = open(mod_pth,'r')
a = f.read()
# run HH-blits
a = a.replace("!HHAPP!", hhblits_app)
a = a.replace("!SEQPTH!", seq_path)
a = a.replace("!LIBDIR!", libdir)
a = a.replace("!OUTA3M!", out_a3m)
a = a.replace("!OUTHHM!", out_hhm)
# align a3m.out to MSA
a = a.replace("!OUTALN!", out_aln)
# convert MSA to a temp a3m file
a = a.replace("!CVTAPP!", convert_aln_app)
a = a.replace("!TMPA3M!", tmp_a3m)
# filter MSA
a = a.replace("!FILTERAPP!", hhfilter_app)
a = a.replace("!FILTEREDA3M!", filtered_a3m)
a = a.replace("!FILTEREDA3MALN!", filtered_aln)
# make pssm
a = a.replace("!PSSMAPP!", pssm_app)
a = a.replace("!PSSMA3M!", pssm_a3m)
a = a.replace("!HHWA3M!", hhw_a3m)
f.close()

g = open(run_prog, 'w+')
g.write(a)
g.close()

os.system("chmod 755 %s" %run_prog)

os.system(queue)
# MYMAX=100
# TOTALMAX=384
# INTERVAL=10
# while True:
#     jobs = int(os.popen('/usr/local/bin/qstat -u sry |wc -l').readline())
#     total_jobs = int(os.popen('/usr/local/bin/qstat |wc -l').readline())
#     if jobs < MYMAX:
#         break
#     elif total_jobs < TOTALMAX:
#         break
#     else:
#         time.sleep(INTERVAL)

os.system("/usr/local/bin/qsub -e %s -o %s -l %s -N %s %s" %(err,out,walltime,tag,run_prog))
print('{run_prog} successfully submitted!'.format(run_prog=run_prog))
time.sleep(1)