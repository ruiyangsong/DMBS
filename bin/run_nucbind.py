#!/usr/bin/python
import os, sys, time, re

pkgdir, datadir, outdir, dataname, pssmpth, hhmpth, sspth, mod_pth, queue = sys.argv[1:]

user = os.popen("whoami").readline().strip()

logdir = '{datadir}/qsub_NucBind_log'.format(datadir=datadir)
os.makedirs(logdir, exist_ok=True)
tag      = "Nuc.{name}".format(name=dataname)
err      = "{logdir}/err".format(logdir=logdir)
out      = "{logdir}/out".format(logdir=logdir)
resource = "walltime=10:00:00"
run_prog = "{logdir}/run_prog.sh".format(logdir=logdir)

## check if the job exists
string = '\n'.join(os.popen('/usr/local/bin/qstat -f | grep -B 1 "Job_Owner.*{user}"'.format(user=user)).readlines())
pattern = re.compile(tag)
result = pattern.findall(string)
if result:
    print("job {tag} already running".format(tag=tag))
    exit(0)


f = open(mod_pth,'r')
a = f.read()
a = a.replace("!PKGDIR!", pkgdir)
a = a.replace("!PSSMPTH!", pssmpth)
a = a.replace("!HHMPTH!", hhmpth)
a = a.replace("!SSPTH!", sspth)
a = a.replace("!OUTDIR!", outdir)
f.close()

g = open(run_prog, 'w+')
g.write(a)
g.close()

os.system("chmod 700 %s" %run_prog)
os.system(queue)

os.system("/usr/local/bin/qsub -e %s -o %s -l %s -N %s %s" %(err,out,resource,tag,run_prog))
print('{run_prog} successfully submitted!'.format(run_prog=run_prog))
time.sleep(1)