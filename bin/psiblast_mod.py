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

# run psi-blast
BLAAPP = "!BLAAPP!"
SEQPTH = "!SEQPTH!"
LIBDIR = "!LIBDIR!"
OUTBLA = "!OUTBLA!"
# align blast.out to MSA
ALNAPP = "!ALNAPP!"
OUTALN = "!OUTALN!"
# convert MSA to a temp blast file
CVTAPP         = "!CVTAPP!"
TMPBLA         = "!TMPBLA!"
# filter MSA
FILTERAPP      = "!FILTERAPP!"
FILTEREDBLA    = "!FILTEREDBLA!"
FILTEREDBLAALN = "!FILTEREDBLAALN!"
# make pssm
PSSMAPP = "!PSSMAPP!"
PSSMBLA = "!PSSMBLA!"
HHWBLA  = "!HHWBLA!"

#
# mian program begins here
#
run_code = 0
## run psi-blast
print('\n@run psi-blast')
run_code += os.system('{app} -query {seq} -db {lib} -out {out_bla} -num_iterations 3'
                      .format(app=BLAAPP, seq=SEQPTH, lib=LIBDIR,out_bla=OUTBLA))
if run_code != 0:
    exit('[ERROR!] run psi-blast failed')

## align blast.out to msa, the WARNING of this app makes run_code != 0
print('\n@align blast.out to msa')
os.system('{app} -i {out_bla} -o {out_aln} -Q {seq} -psi\n'.format(app=ALNAPP, out_bla=OUTBLA, out_aln=OUTALN, seq=SEQPTH))

## convert msa
print('\n@convert and filter msa')
run_code += os.system('{convert_app} {out_aln} {tmp_bla} {filter_app} {filtered_bla} {filtered_bla_aln}'
                      .format(convert_app=CVTAPP, out_aln=OUTALN, tmp_bla=TMPBLA, filter_app=FILTERAPP, filtered_bla=FILTEREDBLA, filtered_bla_aln=FILTEREDBLAALN))
if run_code != 0:
    exit('[ERROR!] convert msa failed')

## make pssm
print('\n@make pssm')
run_code += os.system('{pssm_app} -i {filtered_bla_aln} -P {pssm_bla} -H {hhw_bla}'
                      .format(pssm_app=PSSMAPP, filtered_bla_aln=FILTEREDBLAALN, pssm_bla=PSSMBLA, hhw_bla=HHWBLA))
if run_code != 0:
    exit('[ERROR!] run make pssm failed')

end_time = os.popen('date +"%Y-%m-%d %H:%M:%S"').readline().strip()
print("end at: {end_time}".format(end_time=end_time))