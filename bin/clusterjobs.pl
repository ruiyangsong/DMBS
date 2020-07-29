#!/usr/bin/perl -w
use strict;

my $dname="DMBS";
my $datadir="/public/home/yangserver/DMBS/output/";
my @lst=`cat /public/home/yangserver/DMBS/log/$dname\_running.lst`;
chomp(@lst);
foreach my $l(@lst)
{
    if($l =~ /(\S+)\s+(\S+)/)
    {
	my $id=$1;
	my $target=$2;
	#print "$id,$target\n";
	my $running=`/usr/local/bin/qstat -u yangserver`;
	if (!-e "$datadir/$id")
	{
		`mkdir $datadir/$id`;
	}
	if(-s "$datadir/$id/seq.fasta")
	{
		if(!-s "$datadir/$id/submitted.txt")
	    		{
                	`chmod 777 $datadir/$id/*`;
                    open(ST, ">$datadir/$id/submitted.txt");
                    my $time=`date`;
                	print ST "$target submitted at $time\n";
                	close(ST);
			}
		if($running !~ /SNPs\_$id/)
		{ 
			if(!-s "/public/home/yangserver/$dname/output/$id/output.txt")
		{
			#print "runn\n";
		    `/public/home/yangserver/$dname/bin/run_DMBS.pl $id ` ;
		}
		}
	    }
		if( -s "/public/home/yangserver/$dname/output/$id/output.txt" && !-s "/public/home/yangserver/$dname/output/$id/finished.txt")
		{
			`/public/home/yangserver/$dname/bin/output.pl $id `; 
		}
	}
    }

