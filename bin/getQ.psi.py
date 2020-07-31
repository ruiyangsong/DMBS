#!/usr/bin/python
import os, sys, re, time

JOBTAGPRE=sys.argv[1] # prefix of special Job_Names that you want to monitor

USER="sry"
MYMAX=25     ##change this to the maximum number of jobs by you
TOTALMAX=32  ##change this to the maximum number of jobs by all users
INTERVAL=30
JOBPERNODE=2  ##change this to the maximum number of special jobs by all users on each node

## queue overall jobs
while True:
    jobs = int(os.popen('/usr/local/bin/qstat -u {USER} |wc -l'.format(USER=USER)).readline())
    total_jobs = int(os.popen('/usr/local/bin/qstat |wc -l').readline())
    if jobs <= MYMAX:
        break
    elif total_jobs < TOTALMAX:
        break
    else:
        time.sleep(INTERVAL)

## queue special jobs
while True:
    pattern = re.compile(JOBTAGPRE)
    # up_nodes  = [line.split()[0] for line in os.popen('/usr/local/bin/pbsnodes -l up').readlines()]
    free_nodes = [line.split()[0] for line in os.popen('/usr/local/bin/pbsnodes -l free').readlines()]

    info_str  = os.popen('pbsnodes -a | grep -E "^node| jobs"').readlines()
    for idx in range(len(info_str)):
        if idx % 2 == 0:
            node  = info_str[idx].strip()
            if node in free_nodes:
                ids   = [id.split('/')[-1].strip() for id in info_str[idx+1].split(',')]
                names = [os.popen('/usr/local/bin/qstat -f {id} | grep Job_Name'.format(id=id)).readline().split('=')[-1].strip() for id in ids]
                if len(pattern.findall('.'.join(names))) < JOBPERNODE:
                    print(node)
                    exit(0)
    time.sleep(INTERVAL)

    # hostname  = [info_str[idx].strip() for idx in range(len(info_str)) if idx % 2 == 0]
    # job_ids   = list(map(lambda lst: [id.split('/')[-1].strip() for id in lst], [info_str[idx+1].split(',') for idx in range(len(info_str)) if idx % 2 == 0]))
    # job_names = list(map(lambda lst: [os.popen('/usr/local/bin/qstat -f {id} | grep Job_Name'.format(id=id)).readline().split('=')[-1].strip() for id in lst], job_ids))