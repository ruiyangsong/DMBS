#!/usr/bin/python

import os, sys, time, re

blast_app, libdir, seq_path, outdir, data_name, out_blast, queue,\
alignblast_app, out_aln,\
convert_aln_app, tmp_bla,\
hhfilter_app, filtered_blast, filtered_aln,\
pssm_app, pssm_blast, hhw_blast, \
mod_pth = sys.argv[1:]

logdir = '{outdir}/qsub_blast_log'.format(outdir=outdir)
os.makedirs(logdir, exist_ok=True)
tag_pre  = "psi."
tag      = tag_pre + "{name}".format(name=data_name)
err      = "{logdir}/err".format(logdir=logdir)
out      = "{logdir}/out".format(logdir=logdir)
resource = "walltime=10:00:00"
run_prog = "{logdir}/run_prog.sh".format(logdir=logdir)


## check if the job exists
user = os.popen("whoami").readline().strip()
string = '\n'.join(os.popen('/usr/local/bin/qstat -u {user}'.format(user=user)).readlines())
pattern = re.compile(" "+tag+" ")
result = pattern.findall(string)
if result:
    print("job {tag} already running".format(tag=tag))
    exit(0)


f = open(mod_pth,'r')
a = f.read()
# run psi-blast
a = a.replace("!BLAAPP!", blast_app)
a = a.replace("!SEQPTH!", seq_path)
a = a.replace("!LIBDIR!", libdir)
a = a.replace("!OUTBLA!", out_blast)
# align blast.out to MSA
a = a.replace("!ALNAPP!", alignblast_app)
a = a.replace("!OUTALN!", out_aln)
# convert MSA to a temp blast file
a = a.replace("!CVTAPP!", convert_aln_app)
a = a.replace("!TMPBLA!", tmp_bla)
# filter MSA
a = a.replace("!FILTERAPP!", hhfilter_app)
a = a.replace("!FILTEREDBLA!", filtered_blast)
a = a.replace("!FILTEREDBLAALN!", filtered_aln)
# make pssm
a = a.replace("!PSSMAPP!", pssm_app)
a = a.replace("!PSSMBLA!", pssm_blast)
a = a.replace("!HHWBLA!", hhw_blast)
f.close()

g = open(run_prog, 'w+')
g.write(a)
g.close()

os.system("chmod 755 %s" %run_prog)

hostname = os.popen("{queue} {tag_pre}".format(queue=queue,tag_pre=tag_pre)).readline().strip()
print('available node for psiblast:', hostname)
resource += ",nodes={hostname}".format(hostname=hostname)
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

os.system("/usr/local/bin/qsub -e %s -o %s -l %s -N %s %s" %(err,out,resource,tag,run_prog))
print('{run_prog} successfully submitted!'.format(run_prog=run_prog))
time.sleep(5)