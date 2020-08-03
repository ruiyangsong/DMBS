#!/usr/bin/python
import os, sys, re, time

if len(sys.argv) == 1:
    print("Usage: {filename} [tag_prefix]".format(filename=sys.argv[0]))
    exit(0)

JOBTAGPRE=sys.argv[1] # prefix of special Job_Names that you want to monitor

USER="sry"
MYMAX=80     ##change this to the maximum number of jobs by you
TOTALMAX=384  ##change this to the maximum number of jobs by all users
INTERVAL=30
JOBPERNODE=2  ##change this to the maximum number of special jobs by all users on each node

## Queue overall jobs
while True:
    jobs = int(os.popen('/usr/local/bin/qstat -u {USER} |wc -l'.format(USER=USER)).readline())
    total_jobs = int(os.popen('/usr/local/bin/qstat |wc -l').readline())
    if jobs <= MYMAX:
        break
    elif total_jobs < TOTALMAX:
        break
    else:
        time.sleep(INTERVAL)

## Queue special jobs. They condition is that each node can only submit "JOBPERNODE" special jobs
while True:
    pattern = re.compile(JOBTAGPRE)
    # up_nodes  = [line.split()[0] for line in os.popen('/usr/local/bin/pbsnodes -l up').readlines()]
    free_nodes = [line.split()[0] for line in os.popen('/usr/local/bin/pbsnodes -l free').readlines()]

    if not free_nodes:
        time.sleep(INTERVAL)
        continue

    info_str  = os.popen('/usr/local/bin/pbsnodes -a | grep -E "^node| jobs"').readlines()

    if info_str[-1].strip()[:4] == "node" and info_str[-1].strip() in free_nodes:
        node = info_str[-1].strip()
        print(node)
        exit(0)

    for idx in range(len(info_str) - 1):
        if info_str[idx].strip()[:4] == "node" and info_str[idx + 1].strip()[:4] == "jobs":
            node = info_str[idx].strip()
            if node in free_nodes:
                ids = [id.split('/')[-1].strip() for id in info_str[idx + 1].split(',')]
                names = [os.popen('/usr/local/bin/qstat -f {id} | grep Job_Name'.format(id=id)).readline().split('=')[-1].strip() for id in ids]
                if len(pattern.findall('.'.join(names))) < JOBPERNODE:
                    print(node)
                    exit(0)
        elif info_str[idx].strip()[:4] == "node" and info_str[idx + 1].strip()[:4] == "node":
            node = info_str[idx].strip()
            if node in free_nodes:
                print(node)
                exit(0)
        else:
            continue

    time.sleep(INTERVAL)