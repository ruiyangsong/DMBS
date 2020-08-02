#!/usr/bin/python
import os, sys, re, time

if len(sys.argv) == 1:
    print("Usage: {filename} [pattern_to_grep]".format(filename=sys.argv[0]))
    exit(0)

PATTERN=sys.argv[1] # pattern of the master Job_Names that you want to monitor

USER="sry"
MYMAX=80     ##change this to the maximum number of jobs by you
TOTALMAX=384  ##change this to the maximum number of jobs by all users
INTERVAL=30

MASTERMAX=2  ##change this to the maximum number of master jobs by you

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

## Queue special jobs. They condition is that submitting Master jobs no more than MASTERMAX
while True:
    pattern = re.compile(PATTERN)#pattern.findall('.'.join(names))
    job_cnt = int(os.popen('/usr/local/bin/qstat -u {USER} | grep -E "{PATTERN}" | wc -l'.format(USER=USER,PATTERN=PATTERN)).readline())
    # print(job_cnt)
    if job_cnt < MASTERMAX:
        exit(0)
    else:
        time.sleep(INTERVAL)