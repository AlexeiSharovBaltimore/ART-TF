#!/usr/local/bin/perl

#geneset_match.pl regul_all_targets_genesets_add2022.txt geneset_TRRUST_dir.txt out_ART-TF_TRRUST_dir.txt
#geneset_match.pl regul_all_targets_genesets_add2022.txt signif_Nov2018_2fold_MedEmerCag_genesets.txt out_reg_targets_all_2fold.txt

use strict;
my $usage = "$0 input output\n";
my $arg=0;
my $input_file1 = $ARGV[$arg++] or die $usage;
my $input_file2 = $ARGV[$arg++] or die $usage;
my $output_file = $ARGV[$arg++] or die $usage;

my %geneset;
open (INFO, $input_file1) or die $!;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line=~/^[#\t]/){ next; }
	my($name,$junk,@genes)=split(/\t/, $line);
	$name =~ s/, //;
	$geneset{$name} = \@genes;
}
close INFO;
open (INFO, $input_file2) or die $!;
open (OUT, ">$output_file") or die $!;
my $output_table;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line=~/^[#\t]/){ next; }
	my($name,$junk,@genes)=split(/\t/, $line);
	$name =~ s/, //;
	my $ref = $geneset{$name};
	if(!$ref){ next; }
	my $N1 = @$ref;
	my $N2 = @genes;
	my $Nmin = $N1;
	if($N2<$N1){ $Nmin = $N2; }
	my $N12=0;
	my %hash;
	my @gene_overlap;
	foreach my $gene (@$ref){
		$hash{$gene}=1;
	}
	foreach my $gene (@genes){
		if($hash{$gene}){
			$N12++;
			push(@gene_overlap,$gene);
		}
	}
	$output_table .= "$name\t\t$N12\t".join(",",@gene_overlap)."\n";
	print OUT "$name\t$N1\t$N2\t$N12\n";
}
close INFO;
#print OUT "\n$output_table";
close OUT;
exit(0);


