#!/usr/local/bin/perl

#compare_indirect_direct.pl target_overlap_add_far.txt freq_indirect_direct.txt

use strict;
my $usage = "$0 input_coord TSScoord output\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $!;
my $output_file = $ARGV[$arg++] or die $!;

my $Ngenes = 25090;
my @freqLogPdir;
my @freqLogPindir;
my $log10 = log(10);
open (INFO, $input_file) or die "$input_file NOT found\n";
open (OUT, ">$output_file") or die $!;
while(my $line = <INFO>){
	chop $line;
	my $TF = $line;
	my $change = $line;
	$TF =~ s/_.+$//;
	my $template = $TF."_\\d+";
	my $line = <INFO>;
	chop $line;
	if(!$line){ last; }
	my ($junk,@headers)=split(/\t/, $line);
	my $ndirect = 0;
	while($headers[$ndirect] =~ /$template/){ $ndirect++; }
	my @M;
	my @NNN;
	while(my $line = <INFO>){
		chop $line;
		if(!$line){ last; }
		my ($junk,@data)=split(/\t/, $line);
		my $i1 = @M;
		push(@NNN,$data[$i1]);
		push(@M,\@data);
	}
	if($ndirect==@NNN){ next; }
	my ($Mdir,$Mndir,$n1,$n2)=(0,0,0,0);
	for(my $i=0; $i<$ndirect; $i++){
		if($NNN[$i]<10){ next; }
		for(my $j=$i+1; $j<@NNN; $j++){
			if($NNN[$j]<10){ next; }
			my $N1 = $NNN[$j];
			my $p = $NNN[$i]/$Ngenes;
			if($p==1){ $p=0.99999; }
			if($p==0){ $p=0.00001; }
			my $q = $M[$i]->[$j]/$N1;
			#use hypergeometric distribution test
			my $z = ($q-$p)/sqrt($p*(1-$p)*($Ngenes-$N1)/($Ngenes-1)/$N1);
			if($z<0){ $z=0; }
			#my $p = 1;
			#if($z>0){
			#	$p = 2*(1 - normal_distribution($z));
			#}
			#if($p==0){ print "$z\n"; next; }
			#my $logP = -log($p)/$log10;
			if($j<$ndirect){
				$freqLogPdir[int($z/2)]++;
				$Mdir+=$z; $n1++;
			}else{
				$freqLogPindir[int($z/2)]++;
				$Mndir+=$z; $n2++;
			}
		}
	}
	if($n1>0){ $Mdir/=$n1; }
	if($n2>0){ $Mndir/=$n2; }
	if($n2>0){
		print OUT "$change\t$Mdir\t$Mndir\t$n1\t$n2\n";
	}
}	
close INFO;

for(my $i=0; $i<@freqLogPdir; $i++){
	print OUT "$i\t$freqLogPdir[$i]\t$freqLogPindir[$i]\n";
}
close OUT;
exit(0);

#**********************************************
sub  normal_distribution
#**********************************************
{
my $x = shift;
my $a1=-1.26551223;
my $a2= 1.00002368;
my $a3= 0.37409196;
my $a4= 0.09678418;
my $a5=-0.18628806;
my $a6= 0.27886807;
my $a7=-1.13520398;
my $a8= 1.48851587;
my $a9=-0.82215223;
my $a10=0.17087277;
my $z = abs($x/sqrt(2));
my $t = 1.0 / (1.0 + 0.5 * $z);
my $y = $t*exp(-$z * $z + $a1 + $t * ($a2 + $t * ($a3 + $t * ($a4 + $t * ($a5 + $t *
     ($a6 + $t * ($a7 + $t * ($a8 + $t * ($a9 + $t * $a10)))))))));
if($x < 0.0){ $y = 2.0 - $y; }
$y = 1.0 - 0.5 * $y;
return $y;
}

