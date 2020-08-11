#!/usr/bin/perl -w
use strict;
use List::Util qw/max min/;

my $pkgdir      = "!PKGDIR!";
my $pssmpth     = "!PSSMPTH!";
my $hhmpth      = "!HHMPTH!";
my $sspth       = "!SSPTH!";
my $outdir      = "!OUTDIR!";

if (! -e $outdir){
`mkdir -p $outdir`;
}
chdir "$outdir";

#ss
system("sed '1d' $sspth>1.ss");
my @add=`cat $pkgdir/add`;
open FH1,"1.ss";
open FH2,">2.ss";
print FH2 @add;
while(<FH1>)
{
	chomp;
        s/^\s+//;
        my @l=split/\s+/,$_;
        for(my $i=3; $i<6;$i++){
        print FH2 $l[$i]," ";
        }
        print FH2 "\n";
} 
print FH2 @add;
close FH1;
close FH2;

#pssm
sub std_data
{
         my (@dataset)=@_;
         my $new_elem=0;
         foreach my $elem(@dataset){
         $new_elem=1/(1+exp(-$elem));         
         }
         return ($new_elem);
}

system("sed '1,3d' $pssmpth>pssm-1");
system("tac pssm-1|sed '1,6d'|tac>pssm-2");
my @add2=`cat $pkgdir/add0`; 
open FH1,"pssm-2";
open FH2,">pssm-3";
print FH2 @add2;
while(<FH1>)
{
	chomp;
	s/^\s+//;
	my @l=split/\s+/,$_;
	for(my $i=2; $i<22;$i++){
        my $a=std_data($l[$i]);
	print FH2 $a," ";
	}
	print FH2 "\n";	
} 
print FH2 @add2;
close FH1;
close FH2;

#hhm
open FH1,"$hhmpth";
open FH2,">hhm1";
my $start=0;
while(<FH1>)
{
	chomp;
	if($_ eq "#")
	{
		$start=1;
		next;
	}
	next unless $start;
	s/\*/99999/g;
	print FH2 $_,"\n";
}
close FH1;
close FH2;
system("sed '1,4d' hhm1>hhm2");
system("tac hhm2|sed '1,2d'|tac>hhm3"); 
my @add3=`cat $pkgdir/add1`;
open FH1,"hhm3";
open FH2,">hhm4";
while(<FH1>)
{
	chomp;
	print FH2 substr($_,6),"\t";
	$a=<FH1>;
	print FH2 $a;
	<FH1>;
}
 
close FH1;
close FH2;
open FH1,"hhm4";
open FH2,">1.hhm";
print FH2 @add3;
while(<FH1>)
{
	chomp;
	s/^\s+//;
	my @l=split/\s+/,$_;
	for(my $i=0; $i<20;$i++){
	#print FH2 $l[$i]," ";
	if($l[$i]==99999){print FH2 "0"," ";}
	else{print FH2 2**(-$l[$i]/1000)," ";}
	}
	print FH2 "\n";	
}
print FH2 @add3; 
close FH1;
close FH2;

#combine
my @label=();
my @lpro=`cat ./1.ss`;
foreach my $la(@lpro){
	push(@label,"0");
}
print @label+0,"\n";
open FH1,">spdh";
my @lst1=`cat ./pssm-3`;chomp(@lst1);
my @lst2=`cat ./2.ss`;chomp(@lst2);
my @lst3=`cat ./1.hhm`;chomp(@lst3);
print @lst1+0,"\t",@lst2+0,"\t",@lst3+0,"\n";
for(my $i=0; $i<@lst1; $i++)
{
print FH1 $lst1[$i],"\t",$lst2[$i],"\t",$lst3[$i],"\n";
}
close FH1;

#window
open File,">spdh-pre";
my @pdh=`cat ./spdh`;chomp(@pdh);
print "PDH",@pdh+0,"\n";
my $w=5; 
for(my $i=10; $i<@pdh-10; $i++){ 
	if($label[$i-10]==1||$label[$i-10]==0){
		print File $label[$i-10],"\t",@pdh[$i-5..$i+5],"\n";
	}
 }
close File;	

#forlib
open File,">spdh-pre-lib";
my @spdh=`cat ./spdh-pre`;
foreach my $line(@spdh){
    my @l=split(/\s+/,$line);
    my $i;
    print File $l[$i=0],"\t";
    for($i=1;$i<@l;$i++){
        print File "$i:$l[$i]\t";
    }
    print File "\n";
}
close File;

#SVMpred
system("$pkgdir/svm-scale -r $pkgdir/scale-spdh-MCC  spdh-pre-lib>spdh-pre-5D_scale");
system("$pkgdir/svm-predict  -b 1 spdh-pre-5D_scale $pkgdir/model_spdh_AUC 5D-out");

system("$pkgdir/svm-scale -r $pkgdir/scale-spdh-MCC-R  spdh-pre-lib>spdh-pre-5R_scale");
system("$pkgdir/svm-predict  -b 1 spdh-pre-5R_scale $pkgdir/model_spdh_AUC-R 5R-out");
=pod
#SSITE prediction
my @Bcoach=`cat ./ssite/Bsites_fpt.dat`;chomp(@Bcoach);
my @Pcoach=`cat ./ssite/Bsites_prob_fpt.dat`;chomp(@Pcoach);

die "There is sth wrong \n" if(@Bcoach<1);

my @Pite=split(/\s+/,$Pcoach[0]);
my $length=@Pite;

my @outB=();
my @outP=();
for(my $i=0;$i<$length;$i++){
	push(@outB,"0\n");
	push(@outP,"$Pite[$i]\n");
}

my @llll=split(/\s+/,$Bcoach[0]);
if($llll[5]=~/(\S+)/){
	my $site=$1;
	my @num=split(/\,/,$site);chomp(@num);
	for(my $i=0;$i<@num;$i++){
		my @removed=splice @outB,$num[$i]-1,1,"1\n";
	}
}
chomp(@outB);
chomp(@outP);
open File,">coach-out";
open File2,">Pcoach-out";
for(my $i=0;$i<$length;$i++){
	if($outB[$i]==1){
		if($outP[$i]<0.5){
			print File2 "0.51\n";
		}else{
			print File2 "$outP[$i]\n";
		}
		
	}else{
		print File2 "$outP[$i]\n";
	}
	print File "$outB[$i]\n";
}
close File;
close File2;

#coach-out is about DNA or RNA
my @Ligand=`cat /home/suhong/library/NucModel.lst`;chomp(@Ligand);
my @temp=`cat ./ssite/Bsites_fpt.clr`;chomp(@temp);
my @coachTemplate=();
my $sum_DNA=0;
my $sum_RNA=0;
foreach my $l(@temp){
    if($l=~/(.+)_(NUC)/){
		my $ligand=$1;
		my $num=grep/^$ligand/i,@Ligand;
		if($num>=1){
			foreach my $element(@Ligand){
				my @e=split(/\s+/,$element);
				if($ligand eq $e[0]){
					push(@coachTemplate,"$element\n");
					if(($e[1] eq 'DNA_bind') || ($e[1] eq 'DNA_bind/RNA_bind')){
						$sum_DNA+=1; 
					}
					if(($e[1] eq 'RNA_bind') || ($e[1] eq 'DNA_bind/RNA_bind')){
						$sum_RNA+=1; 
					}					
				}				
			}
		}
    }
	last if($l eq 'cluster 1');
}
my $DNA_per=0;
my $RNA_per=0;
if(@coachTemplate>0){
	$DNA_per=$sum_DNA/(@coachTemplate);
}
if(@coachTemplate>0){
	$RNA_per=$sum_RNA/(@coachTemplate);
}
print $sum_DNA,"\t",$sum_RNA,"\n";

###Dcoach-out, Rcoach-out
my @PCOA=`cat ./Pcoach-out`;chomp(@PCOA);
my @COA=`cat ./coach-out`;chomp(@COA);
my $str="";
open FH,">Dcoach-out-bin";
open FH2,">Dcoach-out-pro";
open File,">Rcoach-out-bin";
open File2,">Rcoach-out-pro";
if($DNA_per==$RNA_per){
	for(my $i=0;$i<@COA;$i++){
		$str=substr($PCOA[$i],0,4);
		print FH "$COA[$i]\n";
		print FH2 "$str\n";
		print File "$COA[$i]\n";
		print File2 "$str\n";		
	}
}else{
	if($DNA_per>$RNA_per){
		for(my $i=0;$i<@COA;$i++){
			$str=substr($PCOA[$i],0,4);
			print FH "$COA[$i]\n";
			print FH2 "$str\n";
			print File "0\n";
			print File2 "0.00\n";
		}		
	}else{
		for(my $i=0;$i<@COA;$i++){
			$str=substr($PCOA[$i],0,4);
			print FH "0\n";
			print FH2 "0.00\n";
			print File "$COA[$i]\n";
			print File2 "$str\n";
		}		
	}
}
close FH;
close FH2;
close File;
close File2;
#SVMpred prediction
my @ProD=`cat ./5D-out`;chomp(@ProD);
my $st  ="";
open FH,">5D-out-pro";
open FH2,">5D-out-bin";
my @tt=split(/\s+/,$ProD[0]);
if($tt[2]==1){
	for(my $i=1;$i<@ProD;$i++){
		my @BIN=split(/\s+/,$ProD[$i]);
		$st=substr($BIN[2],0,4);
		print FH $st,"\n";
		if($BIN[2]>0.50){
		print FH2 1,"\n";
		}else{
		print FH2 0,"\n";
		}
	}
}else{

	for(my $i=1;$i<@ProD;$i++){
		my @BIN=split(/\s+/,$ProD[$i]);
		$st=substr($BIN[1],0,4);
		print FH $st,"\n";
		if($BIN[1]>0.50){
		print FH2 1,"\n";
		}else{
		print FH2 0,"\n";
		}
	}

}
close FH;
close FH2;
my @ProR=`cat ./5R-out`;chomp(@ProR);
my $b	="";
open FH,">5R-out-pro";
open FH2,">5R-out-bin";
my @ttR=split(/\s+/,$ProR[0]);
if($ttR[1]==1){
	for(my $i=1;$i<@ProR;$i++){
		my @BINR=split(/\s+/,$ProR[$i]);
		$b=substr($BINR[1],0,4);
		print FH $b,"\n";
		if($BINR[1]>0.50){
		print FH2 1,"\n";
		}else{
		print FH2 0,"\n";
		}
	}	
}else{
	for(my $i=1;$i<@ProR;$i++){
		my @BINR=split(/\s+/,$ProR[$i]);
		$b=substr($BINR[2],0,4);
		print FH $b,"\n";
		if($BINR[2]>0.50){
		print FH2 1,"\n";
		}else{
		print FH2 0,"\n";
		}
	}
}
close FH;
close FH2;
###out
my @DcoachBin=`cat ./Dcoach-out-bin`;chomp(@DcoachBin);
my @DcoachPro=`cat ./Dcoach-out-pro`;chomp(@DcoachPro);
my @RcoachBin=`cat ./Rcoach-out-bin`;chomp(@RcoachBin);
my @RcoachPro=`cat ./Rcoach-out-pro`;chomp(@RcoachPro);
my @DSvmBin  =`cat ./5D-out-bin`;chomp(@DSvmBin);
my @DSvmPro  =`cat ./5D-out-pro`;chomp(@DSvmPro);
my @RSvmBin  =`cat ./5R-out-bin`;chomp(@RSvmBin);
my @RSvmPro  =`cat ./5R-out-pro`;chomp(@RSvmPro);

my @cscore   =`cat ./ssite/Bsites_fpt.dat`;chomp(@cscore);
my @site     =();
if(@cscore>0){
	@site=split(/\s+/,$cscore[0]);
}
open FH,">D-out";
open FH2,">D-out-pro";
open File,">R-out";
open File2,">R-out-pro";
my $cutoff=0.35;
if($site[0]>$cutoff){
	for(my $i=0;$i<@DSvmBin;$i++){
		print FH "$DcoachBin[$i]\n";
		print FH2 "$DcoachPro[$i]\n";
		print File "$RcoachBin[$i]\n";
		print File2 "$RcoachPro[$i]\n";
	}
}else{
	for(my $i=0;$i<@DSvmBin;$i++){
		print FH "$DSvmBin[$i]\n";
		print FH2 "$DSvmPro[$i]\n";
		print File "$RSvmBin[$i]\n";
		print File2 "$RSvmPro[$i]\n";
	}
}
close FH;
close FH2;
close File;
close File2;
=cut