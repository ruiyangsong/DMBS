#!/usr/bin/perl -w
use strict;

my $id=$ARGV[0];
my $remotehost="222.30.48.146:/home/yangserver/DMBS/output/$id/";
my $outdir="/public/home/yangserver/DMBS/output/$id";

`scp $outdir/output.txt $remotehost`;
`scp $outdir/seq.aln.bla $remotehost`;
`scp $outdir/seq.aln.hhm $remotehost`;
open(FH2,">$outdir/finished.txt");
print FH2 "FINISHED";
close(FH2);