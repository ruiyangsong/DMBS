#!/usr/bin/perl

my $a3m_aln=$ARGV[0];
my $temp_a3m=$ARGV[1];
my $hhfilter_app=$ARGV[2];
my $filtered_a3m=$ARGV[3];
my $filtered_a3m_aln=$ARGV[4];

my @rst=`cat $a3m_aln`;
my $j=0;
open(OUT, ">$temp_a3m");
foreach my $r(@rst)
{
    $j++;
    print OUT ">$j\n";
    print OUT "$r";
}
close(OUT);

`$hhfilter_app -i $temp_a3m -o $filtered_a3m -id 99 -cov 70 -qid 70`;
`egrep -v "^>" $filtered_a3m | sed 's/[a-z]//g' > $filtered_a3m_aln`;
