#!/usr/bin/perl -w
use strict;

######################## S-SITE  Version 20140116 ######################################
# S-SITE is a sequence-based approach to protein-ligand binding site prediction.
# Binding-specific sequence profile-profile alignment is used to recognize homologous
# BioLiP templates. The final ligand-binding sites prediction is done based clustering 
# of multiple templates. Detailed description is available at:
#
# J Yang, A Roy, Y Zhang. Protein-ligand binding site recognition using complementary
#  binding-specific substructure comparison and sequence profile alignment,
#  Bioinformatics, 29:2588-2595 (2013).

# It was developed by Jianyi Yang in the Zhang lab of University of Michigan. 
# If you have any questions, please feel free to contact me at 
# yangji@umich.edu or yang0241@e.ntu.edu.sg
########################################################################################





`hostname`=~/(\S+)/;
my $node=$1;
my $datetime=`date`;
my $path=`pwd`;

print "Hostname: $node\n";
print "Started : $datetime\n";
print "Path    : $path\n";




#=pod      
my $pkgdir      = "!PKGDIR!"; #/home/yangji/I-TASSER3.0
my $libdir      = "!LIBDIR!"; #/home/yangji/ITLIB
my $datadir     = "!DATADIR!";
my $outdir0     = "!OUTDIR!";
my $name        = "!NAME!";
my $seqname     = "!SEQNAME!";
my $tag         = "!TAG!";
my $run         = "!RUN!"; #benchmark or real
my $user        = "!USER!";
my $seqidcut    = "!ID_CUT!";

#=cut
my $bindir      = "$pkgdir/COACH";

my $inputdir    = "$datadir";
#my $seqname     = "seq.txt";
my $outdir      = "$inputdir/ssite";

my $BSLdir      = "$libdir/PSSM";
my $bsseq       = "BSITES.fasta"; ##change to nr version
my $libbsr      = "BSITES.bsr";
my $bsfpt       = "BSITES.fpt"; 

my $jsd_score   = "JSD_coach.dat";

################### parameter set ###########################
#DP parameters
my $wss         = 1.0;   #weight for secondary structure
my $wbs         = 2.0;   #weight for binding site residues
my $wshift      = 0.0;   #weight for shift
my $gapo        = -6.0;  #weight for gap open
my $gape        = -0.5;  #weight for gap extension


my $tan_cut        = 1;   #(1-TANI fingerprint) clustering cutoff
my $tan_cut1       = 0.5; #(1-TANI fingerprint) clustering cutoff
my $cluster_dis    = 8.0; #distance cutoff for bsr clustering
#Template selection parameters
my $nTemplates     = 1;           #how many top templates to select
my $nTemplates1    = 10;          #if not enough templates satisfying cutoff, how many top templates to use
my $q_cut          = 0.5;         #overal quality cutoff for selecting templates
my $zcut           = 0.5;         #z-score cut off to be a binding site residue
################### parameter set ###########################


if($run eq "real")
{
    $seqidcut=1;
}

my $blastdir = "$pkgdir/blast/bin";  # Directory containing  BLAST executables
my $db       = "$libdir/nr/nr";      # Non-redundant sequence database
my $ssdir    = "$pkgdir/PSSpred";
my $blast_parser = "parse_blast.pl";

my $query    = "seq.fasta";  #name to used in my program, do not change
my %Qtemplate=();


if(!-s "$bindir/$blast_parser")
{
    die "$bindir/$blast_parser is missed!\n";
}
if(!-s "$BSLdir/$bsseq")
{
    die "$BSLdir/$bsseq is missed!\n";
}
if(!-s "$BSLdir/$bsfpt")
{
    die "$BSLdir/$bsfpt is missed!\n";
}
if(!-s "$BSLdir/$libbsr")
{
    die "$BSLdir/$libbsr is missed!\n";
}
if(!-s "$inputdir/$seqname")
{
    die "$inputdir/$seqname is missed!\n";
}

`mkdir -p $outdir`;

###############Temporary directory for working ################################
my $workdir  = "/tmp/$user/$tag";
`mkdir -p $workdir`;
`rm -rf $workdir/*`;
chdir "$workdir";


`cp $inputdir/$seqname $query`;
`cp $BSLdir/$bsseq .`;
`cp $BSLdir/$bsfpt .`;
`cp $BSLdir/$libbsr .`;
`cp $bindir/$blast_parser .`;
`cp $bindir//bin/ssite .`;
`cp $bindir/bin/clustering .`;
`cp $bindir/bin/NWalign_simple .`;


##run JSD
if(!-s "$inputdir/$jsd_score")
{
    my $mod=`cat $pkgdir/COACH/JSDmod`;
    $mod    =~ s/\!PKGDIR\!/$pkgdir/mg;
    $mod    =~ s/\!LIBDIR\!/$libdir/mg;
    $mod    =~ s/\!DATADIR\!/$inputdir/mg;
    $mod    =~ s/\!USER\!/$user/mg;
    $mod    =~ s/\!PDBID\!/$name/mg;

    open(FH2,">$inputdir/JSD.pl");
    print FH2 "$mod";
    close(FH2);
    `chmod a+x $inputdir/JSD.pl`;	 
    print "run JSD...\n";
    `$inputdir/JSD.pl`;    
}
if(!-s "$inputdir/$jsd_score")
{
    print "$jsd_score missed\n";
    exit;
}
`cp $inputdir/$jsd_score $workdir`;

if(-s "$outdir/sresult.dat")
{
    print "searching has been done before, copy sresult.dat directly\n";
    `cp $outdir/sresult.dat .`;    
}


############################## start working ###############################
##step 1: run psi-blast and psipred to get the MSA and SS of query
##step 3: search for each fragment against the BSITES Profile library(frag-pocket alignment-->extend to both sids-->assess template with binding residues-->Bscore)
##step 4: select templates with high Bscore, transfer template binding residues to query based on alignment 
##        (one tempate may apprear several times, but with different alignments)


my @AA=qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V );
my @AA1=qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V B Z -);

my %AA2index=
    (
     'A'=>'1', 'R'=>'2', 'N'=>'3', 'D'=>'4', 'C'=>'5', 'Q'=>'6', 'E'=>'7', 'G'=>'8', 'H'=>'9', 'I'=>'10', 'L'=>'11', 'K'=>'12',
     'M'=>'13', 'F'=>'14', 'P'=>'15', 'S'=>'16', 'T'=>'17', 'W'=>'18', 'Y'=>'19', 'V'=>'20', 'B'=>'21', 'Z'=>'22', '-'=>'23');

my %AA2=
    (
     'A'=>'1', 'R'=>'2', 'N'=>'3', 'D'=>'4', 'C'=>'5', 'Q'=>'6', 'E'=>'7', 'G'=>'8', 'H'=>'9', 'I'=>'10', 'L'=>'11', 'K'=>'12',
     'M'=>'13', 'F'=>'14', 'P'=>'15', 'S'=>'16', 'T'=>'17', 'W'=>'18', 'Y'=>'19', 'V'=>'20');
my %SS2num=('C'=>1, 'H'=>2, 'E'=>4);

my %multi=();
my %bslib=();
my %bsligname=();
open(LIB, "$libbsr");
while(my $line = <LIB>)
{
    chomp($line);
    if(substr($line, 11, 1) eq 'P')
    {
        last;
    }
    if($line=~/(\S+)\s+\:(.+)\:(\S+)\:(\d+)\:(\S+)/)
    {
        my $rec =$1;
	my $ligname=$2;
        my $bsno=$3;
        my $bres=$5;
        my $site=$rec . "_" . $bsno;
        $bslib{$site}=$bres;
	push(@{$multi{$rec}}, $site);
	$ligname =~ s/\s+//g;
	$bsligname{$site}=$ligname;
    }
}
close(LIB);


my %libfpt=();
open(FPT, "$bsfpt");
while(my $line = <FPT>)
{
    chomp($line);
    if($line=~/(\S+)\s+(\S+)/)
    {
        my $site      = $1;
        my $fpt       = $2;
        $libfpt{$site} = $fpt;
    }
}
close(FPT);

my @seqtxts=`cat $query`;
my $sequence="";
foreach my $seqtxt(@seqtxts)
{
    if($seqtxt=~/\>/)
    {
	next;
    }
    $seqtxt=~s/\s//mg;
    $seqtxt=~s/\n//mg;
    $sequence=$sequence.$seqtxt;
} 
my %seqQ=();
my $Lch=length $sequence;
open(SEQ,">protein.seq");
printf SEQ ">protein\n";
for(my $i=1;$i<=$Lch;$i++)
{
    $a=substr($sequence, $i-1, 1);
    printf SEQ "$a";
    $seqQ{$i}=$a;   #only for check
    if($i==int($i/60)*60)
    {
        print SEQ "\n";
    }
}   
printf SEQ "\n";
close(SEQ);


if(-s "$outdir/sresult.dat")
{
    goto template_sel;
}



###############################################################
#step 1: run psi-blast and psipred to get the MSA and SS of query
###############################################################

###########  run PSSpred to predict the secondary structure  ########################
if(!-s "$inputdir/seq.ss")
{       
    if(-s "$inputdir/mtx" && -s "$inputdir/pssm.txt" && -s "$inputdir/psitmp.chk" && -s "$inputdir/blast.out")
    {
	`cp $inputdir/mtx .`;
	`cp $inputdir/pssm.txt .`;
	`cp $inputdir/psitmp.chk .`;
	`cp $inputdir/blast.out .`;
    }

    print "$pkgdir/PSSpred/mPSSpred.pl $query $pkgdir $libdir\n";    
    print `$pkgdir/PSSpred/mPSSpred.pl $query $pkgdir $libdir`;
    if(!-s "seq.dat.ss")
    {
	print "PSSpred error!\n";
	exit;
    }
    `mv seq.dat.ss seq.ss`;
    if(!-s "$inputdir/mtx")
    {
	`cp blast.out mtx pssm.txt psitmp.chk $inputdir/`;
    }
    `cp seq.ss seq.dat $inputdir/`;    
}

open(SS, "$inputdir/seq.ss");
my $j=0;
my %sec=();
<SS>;
while(my $line=<SS>)
{
    if($line=~/\s*(\d+)\s+\S+\s+(\S+)/)
    {
	my $j=$1;
	my $ps=$2;
	$sec{$j}=$SS2num{$ps};        
    }
}
close(SS);

%sec = &smooth_ss(\%sec, $Lch);

#run psi-blast
if(-s "$inputdir/blast.out")
{
    printf "copy blast.out directly .....\n";
    `cp $inputdir/blast.out .`;
}
else
{
    printf "running psi-blast .....\n";    
    `$blastdir/blastpgp  -b 1000 -j 3 -h 0.001 -d $db -i protein.seq -C psitmp.chk -Q seq.pssm > blast.out`;
     `cp blast.out $inputdir/`;
}

printf "compute frequency .....\n";
my %close_freq = &wfreq("blast.out");
open(PRO, ">seq.profile");
for(my $i=1; $i<=$Lch; $i++)
{
    my $res=substr($sequence, $i-1, 1);
    printf PRO "$res $sec{$i} ";

    foreach my $A(@AA)
    {
	printf PRO "%.4f ", $close_freq{$i, $A};
    }

    printf PRO "\n";
}
close(PRO);




###############################################################
#step 3: search fragments against the BSITES Profile library
###############################################################
#generate list
my $searchlist="g_template_$seqidcut.lst";
#load receptor sequence
my %recseq=();
open(REC, "$bsseq");
while(my $line = <REC>)
{
    chomp($line);
    if($line =~ /(\S+)\s+(\S+)/)
    {    
        $recseq{$1}=$2;
    }
    elsif($line =~ /^>(\S+)/)
    {
	my $pdb=$1;
	my $seq=<REC>;
	chomp($seq);
	$recseq{$pdb}=$seq;
	
    }
}
close(REC);


my @seqkey = keys %recseq;
if($run eq "benchmark" && $seqidcut <1)
{
    if(-s "$inputdir/$searchlist")
    {
	`cp $inputdir/$searchlist .`;
    }
    else
    {	
	print "prepare template list...\n";
	open(LST, ">$searchlist");
	foreach my $s(@seqkey)
	{
	    my $seq1=$recseq{$s};
	    my $seqid = `./NWalign_simple $seq1 $sequence 3`;	
	    chomp($seqid);
	    if($seqid>$seqidcut)
	    {
		next;
	    }
	    else
	    {
		my @all=@{$multi{$s}};
		foreach my $t(@all)
		{
		    print LST "$t\n";
		}
	    }
	}
	close(LST);
	`cp $searchlist $inputdir/`;
    }
}
else
{
    open(LST, ">$searchlist");
    foreach my $s(@seqkey)
    {
	my @all=@{$multi{$s}};
	foreach my $t(@all)
	{
	    print LST "$t\n";
	}
    }
    close(LST);
}


print "searching BSITES library...\n";
`./ssite seq.profile $searchlist $libdir/PSSM/PSSM $wss $wbs $wshift $gapo $gape`;
`cp sresult.dat $outdir/`;



template_sel:;
###############################################################
#step 4: select templates based on searching result
###############################################################


my %quality=();
my %pbsr=();
my %feature=();

open(FIN, "$jsd_score") or die "No jsdscore\n";
my %JSD=();
while(my $line = <FIN>)
{
    if($line =~ /(\d+)\s+\S+\s+(\S+)/)
    {
	my $res=$1+1;
	my $score=$2;
	$JSD{$res}=$score;
        }
}
close(FIN);



&select_templates("sresult.dat");

my @keys = keys %quality;
foreach my$k(@keys)
{
    #print "$k $feature{$k} $pbsr{$k}\n", ;
}
my $n_selected=@keys;
my $relieve_flag=0;
if(@keys<$nTemplates)
{
    $relieve_flag=1;
    %quality=();
    print "selected number of templates are too small ($n_selected), use top $nTemplates1\n";  
    &select_templates1("sresult.dat");
}
@keys = reverse sort{$quality{$a} <=> $quality{$b}} keys %quality;
$n_selected=@keys;
if($n_selected<1)
{
    die "No templates selected!\n";
}

my @selected=@keys;



my @rst=();
### clustering based on fingerprint similarity

if($relieve_flag==1)
{
    ###compute the pairwise similarity between the fingerprint of the selected template ligand
    my @dis=();
    open(SIM, ">finger_dis.dat");
    for(my $i=0; $i<@keys; $i++)
    {            	    
	my $k=$keys[$i];
	#printf SIM "$k\t";
	my @residues=split(/,/, $pbsr{$k});    
	for(my $j=0; $j<@keys; $j++)
	{	
	    if($i==$j)
	    {
		$dis[$i][$j]=0;
	    }
	    elsif($j<$i)
	    {
		$dis[$i][$j]=$dis[$j][$i];
	    }
	    else
	    {
		my $k1   = $keys[$j];
		my $fpt1 = $libfpt{$k};
		my $fpt2 = $libfpt{$k1};
		
		$dis[$i][$j] = 1.0 - &tannimoto($fpt1, $fpt2);
	    }
	    printf SIM "%.3f ", $dis[$i][$j];
	}
	printf SIM "\n";	
    }
    close(SIM);

    $tan_cut=$tan_cut1;
    @rst=`./clustering finger_dis.dat $n_selected $tan_cut`;    
}
else #no clustering
{
    my $info="cluster 0:";
    for(my $i=0; $i<$n_selected; $i++)
    {
	$info .= "$i ";
    }    
    push(@rst, $info);
}


my %cscore=(); ##total c-score of each cluster
my %clusters=();
my $n_cluster=-1;
my %prob=();
my $bsrncut=0.03*$Lch;
if($bsrncut>10)
{
    $bsrncut=10;
}

&vote(\@rst);
&reRank();
`mv Bsites.dat Bsites_fpt.dat`;
`mv Bsites_prob.dat Bsites_prob_fpt.dat`;
`mv Bsites.clr Bsites_fpt.clr`;



### clustering based on distance between the predicted binding residues 
if(-s "$inputdir/model1.pdb")
{
    `cp $inputdir/model1.pdb .`;
    @rst=&dist_cluster("model1.pdb");

    &vote(\@rst);

    &reRank();
    `mv Bsites.dat Bsites_dis.dat`;
    `mv Bsites_prob.dat Bsites_prob_dis.dat`;
    `mv Bsites.clr Bsites_dis.clr`;
}


`cp Bsites* $outdir/`;

&copy_results($datadir, $outdir0) if($datadir ne $outdir0);

#~~~~~~~~~~~~~~~~~~~~~~~ End procedure  ~~~~~~~~~~~~~~~~~~~~~~~~~#

print "Ending...\n";
my $time=`date`;
printf "ending time: $time";
`sync`;
`sync`;
sleep(1);


`rm -fr $workdir`;
print "Done!\n";

sub reRank
{
    my %cscore=();
    my %pred=();
    my $n=0;
    open(IN, "Bsites.dat");
    while(my $line=<IN>)
    {
	chomp($line);
	if($line =~ /(\S+)\t/)
	{
	    $cscore{$n}=$1;
	    $pred{$n}=$line;
	    $n++;
	}
    }
    close(IN);


    my @keys = sort{$cscore{$b} <=> $cscore{$a}} keys %cscore;

    my %liginf=();
    my %ligcount=();
    my %clr=();
    open(IN, "Bsites.clr") or die "clr is missed\n";
    $n=-1;  
    while(my $line=<IN>)
    {
	#print "$line";
	chomp($line);
	
	if($line =~ /cluster/)
	{
	    $n++;

	    if($n>0)
	    {
		my @keys = sort {$ligcount{$b}<=>$ligcount{$a}} keys %ligcount;
		my $inf="";		
		foreach my $k(@keys)
		{
		    $inf .= "$k $ligcount{$k} ";
		}
		my $m=$n-1;
		$inf =~ s/\s+$//;
		$liginf{$m}=$inf;		
	    }


	    %ligcount=();
	    #print "$n\n";
	}
	else
	{
	    push(@{$clr{$n}}, $line);
	    if($line =~ /(\S+)/)
	    {
		my $site = $1;
		my $lig  = $bsligname{$site};

		if(exists $ligcount{$lig})
		{
		    $ligcount{$lig}++;		    
		}
		else
		{
		    $ligcount{$lig}=1;
		}
	    }
	}
    }
    close(IN);

    if($n>=0)
    {
	my @keys = sort {$ligcount{$b}<=>$ligcount{$a}} keys %ligcount;
	my $inf="";		
	foreach my $k(@keys)
	{
	    $inf .= "$k $ligcount{$k} ";
	}
	my $m=$n;
	$inf =~ s/\s+$//;
	$liginf{$m}=$inf;		
    }

    open(OUT, ">Bsites.dat");
    foreach my $k(@keys)
    {
	print OUT $pred{$k} . "\t" . $liginf{$k} . "\n";
    }
    close(OUT);   


    open(OUT, ">Bsites.clr");
    foreach my $k(@keys)
    {
	my @a=@{$clr{$k}};
	print OUT "cluster $k\n";
	foreach my $j(@a)
	{
	    if($j =~ /(\S+)\t(.+)/)
	    {
		my $site=$1;
		my $inf=$2;
		my $lig=$bsligname{$site};
		print OUT "$site\_$lig\t$inf\n";
	    }
	}
    }
    close(OUT);

 
    my @prob=();
    open(IN, "Bsites_prob.dat");
    while(my $line=<IN>)
    {
	chomp($line);
	if($line =~ /\d+/)
	{
	    push(@prob, $line);
	}
    }
    close(IN);


    open(OUT, ">Bsites_prob.dat");
    foreach my $k(@keys)
    {
	print OUT "$prob[$k]\n";
    }
    close(OUT); 
}




### Ligand Chemical similarity
sub tannimoto
{
    my($f1, $f2)=@_; 
    my @f1bin = split(//, $f1);
    my @f2bin = split(//, $f2);
    my $bins  = $#f1bin;
    $bins     = $#f2bin if($#f2bin < $bins);
    my $NA=0; 
    my $NB=0; 
    my $NAB=0;  
    my $Tan=0;
    for(my $i=0; $i<=$bins; $i++)
    {
	$NA++  if($f1bin[$i] ne '0');
	$NB++  if($f2bin[$i] ne '0');
	$NAB++ if(($f1bin[$i] eq $f2bin[$i]) && ($f1bin[$i] ne '0'))
        }
    my $d = ($NA+$NB-$NAB);
    $Tan = sprintf("%4.3f",($NAB/$d)) if($d > 0) ;
    
    return $Tan;
}


sub select_templates
{
    my ($filename)=@_;   

    open(FIN, $filename) or die "$filename is missed!\n";    
    while(my $line1 = <FIN>)
    {
	chomp($line1);

	if($line1 =~ /(\S+)/)
	{	   	    
	    my @inf=split(/:/, $line1, -1);

	    my $site=$inf[0];
	    #if($bslibd{$site} == 1)
	    #{
		#next;
	    #}
	    my $pdb=substr($site, 0, 5);
	    my $bres=$bslib{$site};    	   

	    my $raw     = $inf[1];
	    my $sia     = $inf[2];
	    my $cov     = $inf[3];
	    my $sib     = $inf[4];
	    my $lcov    = $inf[5];
	    my $predBSR = $inf[6];

	    if($predBSR !~ /\d+/) 
	    {
		#print "ignore $site $predBSR\n";
		next;
	    }

	    $predBSR =~ s/,$//;
	    my @pred=split(/,/, $predBSR);
	    my $cons=0;
	    foreach my $p(@pred)
	    {
		if(!exists $JSD{$p})
		{
		    print "wrong res $p\n";
		    next;
		}
		$cons += $JSD{$p};
	    }
	    if(@pred>0)
	    {
		$cons=$cons/@pred;
	    }


	    my $q = 2.0/(1.0 + exp( -0.5*$raw*$cov - 0.5*$sib*$lcov - 0.2*$cons - 0.00*$sia) )-1.0;

	    #print "$site\t$q\n";
	    
	    $Qtemplate{$site}=$q;
    
	    if($q>$q_cut)
	    {   
		my @pBSR        = split(/,/, $predBSR);
		$quality{$site} = $q;
		$pbsr{$site}    = $predBSR;
		$feature{$site} = sprintf "%.2f %.2f %.2f %.2f %.2f %.2f %.2f", $q, $raw, $sia, $cov, $sib, $lcov, $cons;
	    }
	}
    }
    close(FIN);
}


sub select_templates1
{
    my ($filename)=@_;   

    my @allt = sort {$Qtemplate{$b} <=> $Qtemplate{$a}} keys %Qtemplate;
    my %topT=();
    for(my $i=0; $i<$nTemplates1; $i++)
    {
	$topT{$allt[$i]}=1;
    }


    open(FIN, $filename) or die "$filename is missed!\n";    
    while(my $line1 = <FIN>)
    {
	chomp($line1);

	if($line1 =~ /(\S+)/)
	{	   	    
	    my @inf=split(/:/, $line1, -1);

	    my $site=$inf[0];
	    #if($bslibd{$site} == 1)
	    #{
		#next;
	    #}


	    if(!exists $topT{$site})
	    {
		next;
	    }
	    
	    #print "$site\n";
	   
	    my $pdb=substr($site, 0, 5);
	    my $bres=$bslib{$site};    	   

	    my $raw     = $inf[1];
	    my $sia     = $inf[2];
	    my $cov     = $inf[3];
	    my $sib     = $inf[4];
	    my $lcov    = $inf[5];
	    my $predBSR = $inf[6];

	    if($predBSR !~ /\d+/) 
	    {
		#print "ignore $site $predBSR\n";
		next;
	    }

	    $predBSR =~ s/,$//;
	    my @pred=split(/,/, $predBSR);
	    my $cons=0;
	    foreach my $p(@pred)
	    {
		if(!exists $JSD{$p})
		{
		    print "wrong res $p\n";
		    next;
		}
		$cons += $JSD{$p};
	    }
	    if(@pred>0)
	    {
		$cons=$cons/@pred;
	    }

	    my $q=$Qtemplate{$site};         
	    my @pBSR        = split(/,/, $predBSR);
	    $quality{$site} = $q;
	    $pbsr{$site}    = $predBSR;
	    $feature{$site} = sprintf "%.2f %.2f %.2f %.2f %.2f %.2f %.2f", $q, $raw, $sia, $cov, $sib, $lcov, $cons;	
	}
    }
    close(FIN);
}






sub wfreq
{
    my ($blast)=@_;

    `./$blast_parser -fas blast.out -q protein.seq blast.fasta`;
    my %ALN=();
    my $Pcount=0;
    open(ALN,"blast.fasta") || die "Cant open blast.fasta";
    while(my $line=<ALN>)
    {
        chomp($line);
        if($line =~ /^>(\S+)/)
        {
            my $Pname=$1;
            $Pcount++;
            my $Evalue= $1 if($line =~ /E=(\S+)/);
            $ALN{$Pcount, 0}=$Pname;
            $ALN{$Pcount, 1}=$Evalue;
        }
        else
        {
            $line =~ s/X/-/g;  ###replace X by gap
	    $line=~s/J/-/g;
	    $line=~s/O/-/g;
	    $line=~s/U/-/g;


            $ALN{$Pcount, 2}=$line;
        }
    }
    close(ALN);

    my %freq=();
    if($Pcount > 1)
    {
        $Pcount=1000 if ($Pcount > 1000);
        my %seq_weights = ();
	&load_sequence_weights(\%ALN, \$Pcount, \%AA2index, \%seq_weights);
    

	##compute weighted frequency now
	my $PSEUDOCOUNT= 0.0000001;
	%freq = &frquency(\%ALN, \$Pcount, \%AA2index, \%seq_weights, $PSEUDOCOUNT);
    }
    else
    {
	 my @Qres   = split(//, $ALN{1, 2});
	 for(my $j=0; $j<@Qres; $j++)
	 {
	     foreach my $key (@AA)
	     {
		 $freq{$j+1, $key}=0;
	     }	 
	 }
     }

    return %freq;
}

sub frquency
{
    my ($ALN_ref, $Nseq_ref, $AA_ref, $SW_ref, $pseudocount)=@_;
    my %align   = %$ALN_ref; 
    my $Nseq    = $$Nseq_ref;
    my %weights = %$SW_ref;  
    my %AA2in   = %$AA_ref;
    
    my @Qres   = split(//, $align{1, 2});
    my $Ncol   = $#Qres;
    my %res_count=();

   
    my $Qresno=0;
    my %Qmapping=();
    for(my $j=0; $j<=$#Qres; $j++)
    {
        $res_count{$j}=0;
        if($Qres[$j] ne '-')
        {
            $Qresno++;
            $Qmapping{$Qresno}=$j;
        }
    }

    
    my @ARR=();
    for(my $i=1; $i<=$Nseq; $i++)
    {
        my @res=split(//, $align{$i, 2});
        for(my $j=0; $j<=$#res; $j++)
        {
            $ARR[$i][$j]=$res[$j];
        }
    }
    my $AAcount = keys %AA2;
    my %AA_freq=();
    my %sum_seq_weights=();
    my $k=0;

    for(my $j=0; $j<=$Ncol; $j++)
    {
	if($Qres[$j] eq '-')
	{
	    next;
	}
	$k++;	
        foreach my $key (@AA)
        {
            $AA_freq{$k, $key}=0;
        }
	my $w=0;
        for(my $i=1; $i<=$Nseq; $i++)
        {
	    if(!exists $AA2{$ARR[$i][$j]})
	    {
		next;
	    }
	    $w += $weights{$i};	    
            $AA_freq{$k, $ARR[$i][$j]} += $weights{$i}; ##weighted frequency in clolumn $j
        }
        foreach my $key (@AA)
        {
            $AA_freq{$k, $key} = ($AA_freq{$k, $key}+ $pseudocount)/($w + $AAcount * $pseudocount); ##I change it here	    
        }
    }

    return %AA_freq;
}



# Get Seqeuence Weights based on Henikoff-Henikoff schema
sub load_sequence_weights
{
    my ($ALN_ref, $Nseq_ref, $AA_ref, $seq_weights_ref)=@_;

    my %align  = %$ALN_ref;
    my $Nseq   = $$Nseq_ref;
    my %AA2in  = %$AA_ref;
    
    my $weights = 0;

    my %NRC=();my %RC=();my %seen=();
    for(my $i=1; $i<=$Nseq; $i++)
    {
        my @res=split(//, $align{$i, 2});

        for(my $j=0; $j<=$#res; $j++)
        {
            my $AAN=$AA2in{$res[$j]};
            $RC{$j, $AAN}++;      #number of times $AAN appears at the jth column of the alignment
            if(exists $seen{$j, $AAN})
            {
                next;
            }
            else
            {
                $seen{$j, $AAN}=1;
                $NRC{$j}++;  #total number of AA type in the jth column
            }
        }
    }

    for(my $i=1; $i<=$Nseq; $i++)
    {
        my @res=split(//,$align{$i,2});

        my $ali_L=@res;
        for(my $j=0; $j<=$#res; $j++)
        {
            my $AAN=$AA2in{$res[$j]};
            $$seq_weights_ref{$i} += 1/($NRC{$j} * $RC{$j, $AAN});
        }
        $$seq_weights_ref{$i} = sprintf("%6.4f", $$seq_weights_ref{$i}/$ali_L); ##normalization
        $weights += $$seq_weights_ref{$i};
    }
    return $weights;
}


sub smooth_ss
{
    my ($ssref, $len)=@_;

    my %sec=%$ssref;
    
    #smooth single  --x-- => -----
    for(my $i=3; $i<=$len-2; $i++)
    {
	if($sec{$i}==2 || $sec{$i}==4)
	{
	    my $j=$sec{$i};
	    if($sec{$i-2} != $j)
	    {
		if($sec{$i-1} != $j)
		{
		    if($sec{$i+1} != $j)
		    {
			if($sec{$i+2} != $j)
			{
			    $sec{$i}=1;
			}
		    }
		}
	    }
	}
    }



    #   smooth double 
    #   --xx-- => ------
    for(my $i=1; $i<=$len-5; $i++)
    {
	#helix
	if($sec{$i} != 2)
	{
	    if($sec{$i+1} != 2)
	    {
		if($sec{$i+2} == 2)
		{
		    if($sec{$i+3} == 2)
		    {
			if($sec{$i+4} != 2)
			{
			    if($sec{$i+5} != 2)
			    {
				$sec{$i+2}=1;
				$sec{$i+3}=1;
			    }
			}
		    }
		}
	    }
	}
	
	#beta
	if($sec{$i} != 4)
	{
	    if($sec{$i+1} != 4)
	    {
		if($sec{$i+2} ==4)
		{
		    if($sec{$i+3} == 4)
		    {
			if($sec{$i+4} != 4)
			{
			    if($sec{$i+5} != 4)
			    {
				$sec{$i+2}=1;
				$sec{$i+3}=1;
			    }
			}
		    }
		}
	    }
	}
    }
    

    #smooth connect
    for(my $i=1; $i<=$len-2; $i++)
    {
	if($sec{$i} == 2)
	{
	    if($sec{$i+1} != 2)
	    {
		if($sec{$i+2} == 2)
		{
		    $sec{$i+1}=2;
		}
	    }
	}
	elsif($sec{$i} == 4)
	{
	    if($sec{$i+1} != 4)
	    {
		if($sec{$i+2} == 4)
		{
		    $sec{$i+1}=4;
		}
	    }
	}
    }  

    
    return %sec;
}




sub z_scores
{
     my($score_ref, $z_score_ref)=@_;
     my %r_score  =%$score_ref;

     my @keys = keys %r_score;
     my $Lch = @keys;
     my $z1_a  = 0; my $z1_a2=0; my %zscore=();
     foreach my $i(@keys)
     {
         $z1_a  += $r_score{$i};
         $z1_a2 += $r_score{$i}**2;
     }
     $z1_a/=$Lch;
     $z1_a2/=$Lch;
     my $z1_sd = sqrt($z1_a2 - $z1_a**2);

     foreach my $i(@keys)
     {
         $$z_score_ref{$i}=($r_score{$i}-$z1_a)/$z1_sd;
     }
}



sub dist_cluster
{    
    my ($pdbfile)=@_;
	
    open(QUE, "<$pdbfile") or die "$pdbfile file is missed!\n";
    my %xyz=();
    my $i=0;
    my %resatoms=();
    
    while(my $line = <QUE>)
    {
	if($line =~ /^ATOM/)
	{
	    my $atom_name=substr($line, 12, 4);
	    if($atom_name eq ' CA ')
	    {		    
		my $res_no = substr($line, 22, 5); $res_no =~ s/\s+//g;
		my $check  = $res_no .  $atom_name;
		if(!exists $resatoms{$check}) #ignore alternative location of atoms
		{
		    $resatoms{$check}=1;

		    $xyz{$res_no, 1} = substr($line, 30, 8);
		    $xyz{$res_no, 2} = substr($line, 38, 8);
		    $xyz{$res_no, 3} = substr($line, 46, 8);
		}
	    }
	}
    }
    close(QUE);


    my @dis=();
    ##compute distance matrix
    open(SIM, ">bsr_dis.dat");
    for(my $i=0; $i<@selected; $i++)
    {            	    
	my $k=$selected[$i];
	my @residues=split(/,/, $pbsr{$k});    

	my $x=0; my $y=0; my $z=0;
	foreach my $b(@residues)
	{
	    if(!exists $xyz{$b, 1})
	    {
		print "wrong residue $b\n";
		next;
	    }

	    $x += $xyz{$b, 1};
	    $y += $xyz{$b, 2};
	    $z += $xyz{$b, 3};
	}
	$x = $x/@residues; #center of BSRs
	$y = $y/@residues;
	$z = $z/@residues;    

      
	for(my $j=0; $j<@selected; $j++)
	{	
	    if($i==$j)
	    {
		$dis[$i][$j]=0;
	    }
	    elsif($j<$i)
	    {
		$dis[$i][$j]=$dis[$j][$i];
	    }
	    else
	    {
		my $k1   = $selected[$j];
		my @residues1=split(/,/, $pbsr{$k1});
		
		my $x1=0; my $y1=0; my $z1=0;
		foreach my $b(@residues1)
		{
		    if(!exists $xyz{$b, 1})
		    {
			print "wrong residue $b\n";
			next;
		    }
		    
		    $x1 += $xyz{$b, 1};
		    $y1 += $xyz{$b, 2};
		    $z1 += $xyz{$b, 3};
		}
		$x1 = $x1/@residues1; #center of BSRs
		$y1 = $y1/@residues1;
		$z1 = $z1/@residues1;   
		
		$dis[$i][$j] = &dist($x, $y, $z, $x1, $y1, $z1);
	    }
	    printf SIM "%.3f ", $dis[$i][$j];
	}
	print SIM "\n";
    }
    close(SIM); 

    my $n_selected=@selected;
    my @rst=`./clustering bsr_dis.dat $n_selected $cluster_dis`;     
}

sub dist
{
    my ($x, $y, $z, $x1, $y1, $z1)=@_;
    my $dist=0;
     
    $dist=sqrt( ($x-$x1)**2 + ($y-$y1)**2 + ($z-$z1)**2);
    
    return $dist;
}


sub vote
{
    my ($rst_ref)=@_;

    my @rst=@$rst_ref;



    open(BOUT1, ">Bsites.dat");
    open(BOUT2, ">Bsites_prob.dat");
    open(BOUT3, ">Bsites.clr");
    foreach my $r (@rst)
    {
	chomp($r);
	if($r =~ /^cluster\s+(\d+):(.+)/)
	{
	    $n_cluster++;
	    my $c=$1;
	    my $line=$2;
	    $clusters{$n_cluster}=$line;
	    my @a=split(/\s+/, $line);
	    
	    #make consensus prediction in the cluster @a
	    for(my $k=1; $k<=$Lch; $k++)
	    {
		$prob{$k}=0.0;
	    }
	    
	    my $score=0;
	    my $q0=-1;
	    my $pcount=0;
	    print BOUT3 "cluster $n_cluster\n";
	    foreach my $i(@a)
	    {
		$pcount++;
		
		my $site = $keys[$i];
		my @pred = split(/,/, $pbsr{$site});
		print BOUT3 "$site\t$feature{$site} @pred\n";
		
		if($quality{$site}>$q0)
		{
		    $q0=$quality{$site};
		}
		$score+=$quality{$site};
		
		for(my $k=0; $k<@pred; $k++)
		{
		    $prob{$pred[$k]} += 1;
		}
	    }
	    $score = $score/@a;
	    for(my $k=1; $k<=$Lch; $k++)
	    {
		$prob{$k} /= $pcount if $pcount>0;
		printf BOUT2 "%.2f ", $prob{$k};
	    }
	    print BOUT2 "\n";
	    
	    
	    my @pres=();
	    my %zscore=();
	    &z_scores(\%prob, \%zscore);
	    
	    my @keyr = reverse sort{$zscore{$a} <=> $zscore{$b}} keys %zscore;
	    my $count=0;
	    my $jsd=0;
	    
	    foreach my $k(@keyr)
	    {   	
		if($zscore{$k}>$zcut)
		{		
		    push(@pres, $k);
		    $jsd += $JSD{$k};
		    $count++;
		    if($count>=$bsrncut ) #for hard target, make at most min(10, 0.03*$Lch) predictions
		    {
			#last;
		    }
		} 
	    }
	    if($count==0)
	    {
		$zcut=0;
		foreach my $k(@keyr)
		{   	
		    if($zscore{$k}>$zcut)
		    {		
			push(@pres, $k);
			$jsd += $JSD{$k};
			$count++;
			if($count>=3) #for hard target, make at most min(10, 0.03*$Lch) predictions
			{
			    last;
			}	 
		    }    
		}
	    }
	    
	    
	    $jsd = $jsd/$count if $count>0;
	    
	    

	    my $size  = @a/@keys;
	    my $size1 = log(1+@a ** (1.0/3));
	    #$cscore{$n_cluster} = 2/(1+exp( -$size*$score - 0.2*$size1 -0.2*$jsd)  ) -1;
	    $cscore{$n_cluster} = 2/(1+exp( -$size*$q0 - 0.1*$size1 -0.2*$jsd)  ) -1;
	    
	    #$cscore{$k} = 2/(1+exp( -$size*$highestq[$k] - 0.2*$size1 - 0.2*$conser[$k] )  ) -1;
	    
	    @keyr= sort {$a <=> $b} @pres;
	    my $ns=@a;
	    printf BOUT1 "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t", $cscore{$n_cluster}, $q0, $size, $size1, $jsd;
	    print BOUT1 join(",", @keyr) . "\n";
	}
    }
    close(BOUT1);
    close(BOUT2);
    close(BOUT3);    
}

`mv -f $inputdir/JSD* $inputdir/blast.out $inputdir/mtx $inputdir/psitmp.chk $inputdir/pssm.txt $inputdir/seq.dat $inputdir/seq.ss $outdir`;
sub copy_results
{
    my ($dir1, $dir2)=@_;
  
    `mkdir -p $dir2/ssite`;
    chdir "$dir1/ssite";
    `cp -f Bsites_* $dir2/ssite/`;
}
