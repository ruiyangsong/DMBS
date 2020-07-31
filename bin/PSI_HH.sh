#!/bin/bash

[ $# != 1 ] && {
    echo "USAGE: $0 [dataname]"
    exit
}

USER=`whoami`
echo "begin at: `date +"%Y-%m-%d %H:%M:%S"`"
echo "user: $USER"
echo "hostname: `hostname`"
echo -e "pwd: `pwd`\n"

# ----------------------------------------------------------------------------------------------------------------------
# some static directories, input file (seq.fasta, mut_pos.txt and mut_res.txt) should be placed in $OUTDIR in advance
# ----------------------------------------------------------------------------------------------------------------------
INTERVAL=300
COUNTMAX=120
DATANAME=$1
BINDIR=/public/home/$USER/DMBS/src/bin
MDLDIR=/public/home/$USER/DMBS/model
OUTDIR=/public/home/$USER/DMBS/output/$DATANAME
#TMPDIR=/tmp/DMBS/$DATANAME
QUEUE=$BINDIR/getQ.py
QUEUE_PSI=$BINDIR/getQ.psi.py
# psi-blast
PSIBLAST_APP=/public/home/caobx/local/app/blast/bin/psiblast
PSIBLAST_LIB=/library/uniref100/uniref100
# hh-search
HHBLITS_APP=/public/home/caobx/local/app/hhsuite3/bin/hhblits
HHBLITS_LIB=/library/uniclust30/uniclust30_2018_08
HHFILTER_APP=/public/home/caobx/local/app/hhsuite3/bin/hhfilter

# ----------------------------------------------------------------------------------------------------------------------
# main program begins here, check file
# ----------------------------------------------------------------------------------------------------------------------
echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> check seq.fasta, mut_pos.txt and mut_res.txt file at $OUTDIR"
[ -e $OUTDIR/seq.fasta ] || {
    echo '[ERROR!] found no seq.fasta at $OUTDIR'
    exit
}
[ -e $OUTDIR/mut_pos.txt ] || {
    echo '[ERROR!] found no mut_pos.txt at $OUTDIR'
    exit
}
[ -e $OUTDIR/mut_res.txt ] || {
    echo '[ERROR!] found no mut_res.txt at $OUTDIR'
    exit
}

# ----------------------------------------------------------------------------------------------------------------------
# HH-BLITS for native sequence
# ----------------------------------------------------------------------------------------------------------------------
echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> qsub HH-BLITS for native sequence, see the log dir at $OUTDIR for details"
OUTA3M=$OUTDIR/seq.a3m
OUTHHM=$OUTDIR/seq.hhm

OUTA3MALN=$OUTDIR/seq.aln.hhm

TEMPA3M=$OUTDIR/seq.a3m.tmp
FILTEREDA3M=$OUTDIR/seq.filtered.a3m.tmp
FILTEREDA3MALN=$OUTDIR/seq_hhm_filtered.aln

PSSM_A3M=$OUTDIR/seq_hhm.pssm
HHW_A3M=$OUTDIR/seq_hhm.hhw
if [ ! -e $OUTDIR/seq_hhm.pssm ] || [ ! -e $OUTDIR/seq_hhm.hhw ]
then
    $BINDIR/run_hhblits.py $HHBLITS_APP $HHBLITS_LIB $OUTDIR/seq.fasta $OUTDIR $DATANAME $OUTA3M $OUTHHM $QUEUE $OUTA3MALN $BINDIR/convert_aln.pl $TEMPA3M $HHFILTER_APP $FILTEREDA3M $FILTEREDA3MALN $BINDIR/make_pssm.py $PSSM_A3M $HHW_A3M $BINDIR/hhblits_mod.py || {
        echo '[ERROR!] qsub HHBLITS for native sequence failed, check log dir at $OUTDIR for details'
        exit
    }
else
    echo "[WARNING!] qsub HH-BLITS canceled as seq_hhm.pssm and seq_hhm.hhw already exist at $OUTDIR"
fi

# -----------------------------------------------------------------------------------------------------------------------s
# make mutation
# ----------------------------------------------------------------------------------------------------------------------
echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> get POS and MUT from mut_pos.txt and mut_res.txt"
POS=`cat $OUTDIR/mut_pos.txt`
MUT=`cat $OUTDIR/mut_res.txt`

echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> split native fasta file to HEAD and SEQ, and get wild residue type at pos_$POS"
SEQHEADTMP=`head -n 1 $OUTDIR/seq.fasta`
SEQHEAD=`echo $SEQHEADTMP |sed 's/\r//g'`
SEQTMP=`grep -v ">" $OUTDIR/seq.fasta`
SEQ=`echo $SEQTMP |sed 's/\s//g'`
WT=`echo ${SEQ:$POS-1:1}`

# ----------------------------------------------------------------------------------------------------------------------
#  HH-BLITS and PSSM calculation for each mutated sequence
# ----------------------------------------------------------------------------------------------------------------------
PREPARED=""
for MT in `echo $MUT | awk -F "" '{for (i=1;i<=NF;i++){print($i)}}'`
do
    [ $MT = $WT ] && continue
    echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> make mutated sequence for $WT at pos_$POS to $MT"
    PREPARED=$PREPARED$MT
    MUTDIR=$OUTDIR/$MT
    MUTNAME=$DATANAME.$POS.$MT
    [ -e $MUTDIR ] || mkdir -p $MUTDIR
    SEQPRE=`echo ${SEQ:0:$POS-1}`
    SEQSUF=`echo ${SEQ:$POS}`
    SEQMUT=$SEQPRE$MT$SEQSUF
    echo $SEQHEAD > $MUTDIR/mut.fasta
    echo $SEQMUT >> $MUTDIR/mut.fasta

    #
    # HH-BLITS
    #
    echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> qsub HH-BLITS for $MUTNAME sequence, see the log dir at $MUTDIR for details"
    OUTA3M=$MUTDIR/mut.a3m
    OUTHHM=$MUTDIR/mut.hhm

    OUTA3MALN=$MUTDIR/mut.aln.hhm

    TEMPA3M=$MUTDIR/mut.a3m.tmp
    FILTEREDA3M=$MUTDIR/mut.filtered.a3m.tmp
    FILTEREDA3MALN=$MUTDIR/mut_hhm_filtered.aln

    PSSM_A3M=$MUTDIR/mut_hhm.pssm
    HHW_A3M=$MUTDIR/mut_hhm.hhw
    if [ ! -e $MUTDIR/mut_hhm.pssm ] || [ ! -e $MUTDIR/mut_hhm.hhw ]
    then
        $BINDIR/run_hhblits.py $HHBLITS_APP $HHBLITS_LIB $MUTDIR/mut.fasta $MUTDIR $MUTNAME $OUTA3M $OUTHHM $QUEUE $OUTA3MALN $BINDIR/convert_aln.pl $TEMPA3M $HHFILTER_APP $FILTEREDA3M $FILTEREDA3MALN $BINDIR/make_pssm.py $PSSM_A3M $HHW_A3M $BINDIR/hhblits_mod.py || {
            echo '[ERROR!] qsub HH-BLITS for $MUTNAME sequence failed, check log dir at $MUTDIR for details'
            exit
        }
    else
        echo "[WARNING!] qsub HH-BLITS canceled as mut_hhm.pssm and mut_hhm.hhw already exist at $MUTDIR"
    fi
done

# ----------------------------------------------------------------------------------------------------------------------
# PSI-BLAST for native sequence
# ----------------------------------------------------------------------------------------------------------------------
echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> qsub PSI-BLAST for native sequence, see the log dir at $OUTDIR for details"
OUTBLA=$OUTDIR/seq.blast

OUTBLAALN=$OUTDIR/seq.aln.bla

TEMPBLA=$OUTDIR/seq.bla.tmp
FILTEREDBLA=$OUTDIR/seq.filtered.bla.tmp
FILTEREDBLAALN=$OUTDIR/seq_bla_filtered.aln

PSSM_BLA=$OUTDIR/seq_bla.pssm
HHW_BLA=$OUTDIR/seq_bla.hhw
if [ ! -e $OUTDIR/seq_bla.pssm ] || [ ! -e $OUTDIR/seq_bla.hhw ]
then
    $BINDIR/run_psiblast.py $PSIBLAST_APP $PSIBLAST_LIB $OUTDIR/seq.fasta $OUTDIR $DATANAME $OUTBLA $QUEUE_PSI $BINDIR/alignblast.pl $OUTBLAALN $BINDIR/convert_aln.pl $TEMPBLA $HHFILTER_APP $FILTEREDBLA $FILTEREDBLAALN $BINDIR/make_pssm.py $PSSM_BLA $HHW_BLA $BINDIR/psiblast_mod.py || {
        echo '[ERROR!] qsub PSI-BLAST for native sequence failed, check log dir at $OUTDIR for details'
        exit
    }
else
    echo "[WARNING!] qsub PSI-BLAST canceled as seq_bla.pssm and seq_bla.hhw already exist at $OUTDIR"
fi

# ----------------------------------------------------------------------------------------------------------------------
# PSI-BLAST and PSSM calculation for each mutated sequence
# ----------------------------------------------------------------------------------------------------------------------
PREPARED=""
for MT in `echo $MUT | awk -F "" '{for (i=1;i<=NF;i++){print($i)}}'`
do
    [ $MT = $WT ] && continue
    echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> make mutated sequence for $WT at pos_$POS to $MT"
    PREPARED=$PREPARED$MT
    MUTDIR=$OUTDIR/$MT
    MUTNAME=$DATANAME.$POS.$MT
    [ -e $MUTDIR ] || mkdir -p $MUTDIR
    SEQPRE=`echo ${SEQ:0:$POS-1}`
    SEQSUF=`echo ${SEQ:$POS}`
    SEQMUT=$SEQPRE$MT$SEQSUF
    echo $SEQHEAD > $MUTDIR/mut.fasta
    echo $SEQMUT >> $MUTDIR/mut.fasta

    #
    # PSI-BLAST
    #
    echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> qsub PSI-BLAST for $MUTNAME sequence, see the log dir at $MUTDIR for details"
    OUTBLA=$MUTDIR/mut.blast

    OUTBLAALN=$MUTDIR/mut.aln.bla

    TEMPBLA=$MUTDIR/mut.bla.tmp
    FILTEREDBLA=$MUTDIR/mut.filtered.bla.tmp
    FILTEREDBLAALN=$MUTDIR/mut_bla_filtered.aln

    PSSM_BLA=$MUTDIR/mut_bla.pssm
    HHW_BLA=$MUTDIR/mut_bla.hhw
    if [ ! -e $MUTDIR/mut_bla.pssm ] || [ ! -e $MUTDIR/mut_bla.hhw ]
    then
        $BINDIR/run_psiblast.py $PSIBLAST_APP $PSIBLAST_LIB $MUTDIR/mut.fasta $MUTDIR $MUTNAME $OUTBLA $QUEUE_PSI $BINDIR/alignblast.pl $OUTBLAALN $BINDIR/convert_aln.pl $TEMPBLA $HHFILTER_APP $FILTEREDBLA $FILTEREDBLAALN $BINDIR/make_pssm.py $PSSM_BLA $HHW_BLA $BINDIR/psiblast_mod.py || {
            echo '[ERROR!] qsub PSI-BLAST for $MUTNAME sequence failed, check log dir at $MUTDIR for details'
            exit
        }
    else
        echo "[WARNING!] qsub PSI-BLAST canceled as mut_bla.pssm and mut_bla.hhw already exist at $MUTDIR"
    fi
done








## ----------------------------------------------------------------------------------------------------------------------
## This loop for integrating features that pre-calculated by PSI-BLAST and HH-BLITS for each mutation
## ----------------------------------------------------------------------------------------------------------------------
#for MT in `echo $MUT | awk -F "" '{for (i=1;i<=NF;i++){print($i)}}'`;
#do
#    [ $MT = $WT ] && continue
#    MUTDIR=$OUTDIR/$MT
#    TIMEOUT=1
#    while [ $TIMEOUT -le $COUNTMAX ] && ([ ! -e $OUTDIR/seq_bla.pssm ] || [ ! -e $OUTDIR/seq_bla.hhw ] || [ ! -e $OUTDIR/seq_hhm.pssm ] || [ ! -e $OUTDIR/seq_hhm.hhw ] || [ ! -e $MUTDIR/mut_bla.pssm ] || [ ! -e $MUTDIR/mut_bla.hhw ]  || [ ! -e $MUTDIR/mut_hhm.pssm ] || [ ! -e $MUTDIR/mut_hhm.hhw ])
#    do
#        echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> check whether PSI-BLAST and HH-BLITS are done properly, current TIMEOUT is $TIMEOUT"
#        sleep $INTERVAL
#        let TIMEOUT+=1
#    done
#    [ $TIMEOUT -le $COUNTMAX ] || {
#        echo "[TIMEOUT!] wait PSI-BLAST and HH-BLITS timeout, INTERVAL is $INTERVAL, COUNTMAX is $COUNTMAX"
#        exit
#    }
#    echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> integrate pre-calculated sequence-based features"
#    [ -e $MUTDIR/DMBS_feature.csv ] || {
#        $BINDIR/integrate.py $POS $MT $OUTDIR $BINDIR/BLOSUM62.txt || {
#            echo "[ERROR!] integrate features failed, current mutation is $WT at pos_$POS to $MT"
#            exit
#        }
#    }
#done
#
## ----------------------------------------------------------------------------------------------------------------------
## Make prediction
## ----------------------------------------------------------------------------------------------------------------------
#echo "[`date +"%Y-%m-%d %H:%M:%S"`] --> make prediction by model at $MDLDIR"
#[ -e $OUTDIR/pred_all.txt ] || $BINDIR/pred.py $OUTDIR $PREPARED $MDLDIR/var.pickle








echo -e "\nend at: `date +"%Y-%m-%d %H:%M:%S"`"