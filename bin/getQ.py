#!/usr/bin/python
import os
import time
USER="sry"
MYMAX=80     ##change this to the maximum number of jobs by you
TOTALMAX=384  ##change this to the maximum number of jobs by all users
INTERVAL=30
while True:
    jobs = int(os.popen('/usr/local/bin/qstat -u {USER} |wc -l'.format(USER=USER)).readline())
    total_jobs = int(os.popen('/usr/local/bin/qstat |wc -l').readline())
    if jobs <= MYMAX:
        break
    elif total_jobs < TOTALMAX:
        break
    else:
        time.sleep(INTERVAL)