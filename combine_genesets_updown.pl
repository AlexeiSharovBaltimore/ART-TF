#!/usr/local/bin/perl

#combine_genesets_updown.pl regul_all_targets_genesets_add.txt regul_all_targets_genesets_updown.txt

use strict;
my $usage = "$0 inputfile\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $usage;
my $output_file = $ARGV[$arg++] or die $usage;

open (INFO, $input_file) or die $!;
open (OUT, ">$output_file") or die $!;
my %hashSymb;
my $TF_old;
while(my $line = <INFO>){
	chop $line;
	if(!$line || $line =~ /^!|^</){ next; }
	my ($name,$N,@symbols)=split(/\t/, $line);
	if(!$name){ next; }
	my ($TF,$dir) = split(/_/,$name);
	if($TF_old && $TF_old ne $TF){
		my @symbols = keys %hashSymb;
		my $N1 = @symbols;
		print OUT "$TF_old\t$N1\t".join("\t",@symbols)."\n";
		%hashSymb = ();
	}
	for(my $i=0; $i<@symbols; $i++){
		$hashSymb{$symbols[$i]}=1;
	}
	$TF_old = $TF;
}
if($TF_old){
	my @symbols = keys %hashSymb;
	my $N1 = @symbols;
	print OUT "$TF_old\t$N1\t".join("\t",@symbols)."\n";
}
close INPUT;
close OUT;
exit(0);

