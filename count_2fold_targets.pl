#!/usr/local/bin/perl

#count_2fold_targets.pl regul_all_targets_genesets_add.txt signif_Nov2017_2fold_MedEmerCag_genesets.txt regul_all_targets_2fold.txt

use strict;
my $usage = "$0 inputfile\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $usage;
my $signif_genesets = $ARGV[$arg++] or die $usage;
my $output_file = $ARGV[$arg++] or die $usage;

my %hashOriginal=("AHSA2P","AHSA2","ATP5F1A","ATP5A1","ATP5F1AP1","ATP5A1P1","ATP5F1AP2","ATP5A1P2","ATP5F1AP3","ATP5A1P3","ATP5F1AP7","ATP5A1P7","ATP5F1AP8","ATP5A1P8","ATP5F1B","ATP5B","ATP5F1C","ATP5C1","ATP5F1CP1","ATP5C1P1","ATP5F1D","ATP5D","ATP5F1E","ATP5E","ATP5F1EP2","ATP5EP2","ATP5PB","ATP5F1","ATP5PBP2","ATP5F1P2","ATP5PBP3","ATP5F1P3","ATP5PBP5","ATP5F1P5","ATP5PBP6","ATP5F1P6","ATP5PBP7","ATP5F1P7","ATP5MC1","ATP5G1","ATP5MC1P4","ATP5G1P4","ATP5MC1P5","ATP5G1P5","ATP5MC1P7","ATP5G1P7","ATP5MC2","ATP5G2","ATP5MC2P1","ATP5G2P1","ATP5MC2P3","ATP5G2P3","ATP5MC2P4","ATP5G2P4","ATP5MC3","ATP5G3","ATP5PD","ATP5H","ATP5PDP2","ATP5HP2","ATP5PDP3","ATP5HP3","ATP5PDP4","ATP5HP4","ATP5ME","ATP5I","ATP5PF","ATP5J","ATP5MF","ATP5J2","ATP5MFP4","ATP5J2P4","ATP5MFP5","ATP5J2P5","ATP5MF-PTCD1","ATP5J2-PTCD1","ATP5PFP1","ATP5JP1","ATP5MG","ATP5L","ATP5MGL","ATP5L2","ATP5MGP2","ATP5LP2","ATP5MGP3","ATP5LP3","ATP5MGP4","ATP5LP4","ATP5MGP5","ATP5LP5","ATP5PO","ATP5O","ATP5IF1","ATPIF1","DEPP1","C10orf10","ATP5MPL","C14orf2","CFAP97D1","C17orf105","RMC1","C18orf8","MIR29B2CHG","C1orf132","RHEX","C1orf186","ODR4","C1orf27","MIR1-1HG-AS1","C20orf166-AS1","SNORC","C2orf82","CAND1","CAND1.11","EFL1P1","EFTUD1P1","EP400P1","EP400NL","GUCY1A1","GUCY1A3","GUCY1B1","GUCY1B3","TRMT9B","KIAA1456","IQCN","KIAA1683","SMIM37","LINC00116","MIR570HG","LINC00969","ROCR","LINC02095","MIR4453HG","LINC02486","RNF227","LINC02581","SLC49A3","MFSD7","NUDT4B","NUDT4P1","PRSS46P","PRSS46","COP1","RFWD2","MTREX","SKIV2L2","TBXT","T","STIMATE","TMEM110","CLTRN","TMEM27","RXYLT1","TMEM5","MACO1","TMEM57","ATP5MD","USMG5","ATP5MDP1","USMG5P1","RAB6D","WTH3DI","ZRSR2P1","ZRSR1");

my %hashSignif;
open (INFO, $signif_genesets) or die $!;
while(my $line = <INFO>){
	chop $line;
	if(!$line || $line =~ /^!|^</){ next; }
	my ($name,$N,@symbols)=split(/\t/, $line);
	if(!$name){ next; }
	for(my $i=0; $i<@symbols; $i++){
		if($hashOriginal{$symbols[$i]}){ $symbols[$i]=$hashOriginal{$symbols[$i]}; }
	}
	$hashSignif{$name} = \@symbols;
}
close INFO;

open (INFO, $input_file) or die $!;
open (OUT, ">$output_file") or die $!;
print OUT "ChIP-seq\tInduced_TF\n";
while(my $line = <INFO>){
	chop $line;
	if(!$line || $line =~ /^!|^</){ next; }
	my ($name,$N,@symbols)=split(/\t/, $line);
	if(!$name){ next; }
	#while($line = <INFO>){
	#	if($line =~ /^\tlogratio\t/){ last; }
	#}
	#chop $line;
	#my ($junk1,$type,@logratio)=split(/\t/, $line);
	my $N2fold=0;
	my $ref = $hashSignif{$name};
	if($ref){
		my %hash;
		foreach my $symb (@$ref){ $hash{$symb}=1; }
		foreach my $symb (@symbols){
			if ($hash{$symb}){ $N2fold++; }
		}
	}
	#for(my $i=0; $i<@logratio; $i++){
	#	if($logratio[$i]<-0.301 || $logratio[$i]>0.301){ $N2fold++; }
	#}
	print OUT "$name\t$N\t $N2fold\n";
}
close INPUT;
close OUT;
exit(0);

