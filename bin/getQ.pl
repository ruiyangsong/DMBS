#!/usr/bin/perl -w
use strict;

my $user=`whoami`;
$user =~ s/\n//;
my $MYMAX=100;       ###change this to the maximum number of jobs by you
my $TOTMAX=384;      ##change this to the maximum number of jobs by all users
my $INTERVAL=10;

while(1)
{
    my $count="";
    while($count eq "")
    {
	$count=`qstat -u $user |wc -l`;
	$count =~ s/\s+//g;
    }
    
    my $q="default";   
    

    if($count<$MYMAX)
    {
	#print "$q";
	exit(0);
    }
    else
    {
	$count="";
	while($count eq "")
	{
	    $count=`qstat |wc -l`;
	    $count =~ s/\s+//g;
	}

	if($count<$TOTMAX)
	{
	    #print "$q";
	    exit(0);
	}
	else
	{
	    sleep($INTERVAL);
	}
    }
}
