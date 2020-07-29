#!/usr/bin/python
import os

begin_time = os.popen('date +"%Y-%m-%d %H:%M:%S"').readline().strip()
user       = os.popen("whoami").readline().strip()
hostname   = os.popen('hostname').readline().strip()
pwd        = os.popen('pwd').readline().strip()
print("begin at: {begin_time}".format(begin_time=begin_time))
print("user: {user}".format(user=user))
print("hostname: {hostname}".format(hostname=hostname))
print("pwd: {pwd}".format(pwd=pwd))

# run HH-blits
HHAPP  = "!HHAPP!"
SEQPTH = "!SEQPTH!"
LIBDIR = "!LIBDIR!"
OUTA3M = "!OUTA3M!"
OUTHHM = "!OUTHHM!"
# align a3m.out to MSA
OUTALN = "!OUTALN!"
# convert MSA to a temp a3m file
CVTAPP = "!CVTAPP!"
TMPA3M = "!TMPA3M!"
# filter MSA
FILTERAPP      = "!FILTERAPP!"
FILTEREDA3M    = "!FILTEREDA3M!"
FILTEREDA3MALN = "!FILTEREDA3MALN!"
# make pssm
PSSMAPP = "!PSSMAPP!"
PSSMA3M = "!PSSMA3M!"
HHWA3M  = "!HHWA3M!"

#
# mian program begins here
#
run_code = 0
## run HH-blits
print('\n@run HH-blits')
run_code += os.system('{app} -i {seq} -d {lib} -oa3m {out_a3m} -ohhm {out_hhm} -n 3 -maxfilt 500000 -diff inf -id 99 -cov 70'
                      .format(app=HHAPP, seq=SEQPTH, lib=LIBDIR,out_a3m=OUTA3M, out_hhm=OUTHHM))
if run_code != 0:
    exit('[ERROR!] run HH-blits failed')

## align a3m to msa
print('\n@align a3m to MSA')
run_code += os.system('egrep -v "^>" {out_a3m} | sed \'s/[a-z]//g\' > {out_aln}'.format(out_a3m=OUTA3M,out_aln=OUTALN))
if run_code != 0:
    exit('[ERROR!] align a3m to msa failed')

## convert msa
print('\n@convert and filter msa')
run_code += os.system('{convert_app} {out_aln} {tmp_a3m} {filter_app} {filtered_a3m} {filtered_a3m_aln}'
                      .format(convert_app=CVTAPP, out_aln=OUTALN, tmp_a3m=TMPA3M, filter_app=FILTERAPP, filtered_a3m=FILTEREDA3M, filtered_a3m_aln=FILTEREDA3MALN))
if run_code != 0:
    exit('[ERROR!] convert msa failed')

## make pssm
print('\n@make pssm')
run_code += os.system('{pssm_app} -i {filtered_a3m_aln} -P {pssm_a3m} -H {hhw_a3m}'
                      .format(pssm_app=PSSMAPP, filtered_a3m_aln=FILTEREDA3MALN, pssm_a3m=PSSMA3M, hhw_a3m=HHWA3M))
if run_code != 0:
    exit('[ERROR!] make pssm failed')

end_time = os.popen('date +"%Y-%m-%d %H:%M:%S"').readline().strip()
print("end at: {end_time}".format(end_time=end_time))