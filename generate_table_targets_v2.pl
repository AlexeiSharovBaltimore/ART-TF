#!/usr/local/bin/perl

#generate_table_targets_v2.pl regul_dir_targets_genesets_add.txt regul_all_targets_genesets_add.txt output.txt

use strict;
my $usage = "$0 input output\n";
my $arg=0;
my $input_file1 = $ARGV[$arg++] or die $usage;
my $input_file2 = $ARGV[$arg++] or die $usage;
my $output_file = $ARGV[$arg++] or die $usage;

my %direct;
my %geneset;
my $name_old;
my @data;
open (INFO, $input_file1) or die $!;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line=~/^#/){ next; }
	my($name,$junk,@genes)=split(/\t/, $line);
	my ($TF,$dir) = split(/_/,$name);
	$direct{$TF}++;
	if($name && $name_old){
		$geneset{$name_old} = [@data];
		@data = ();
	}
	push(@data,\@genes);
	if($name){
		$name_old = $name;
	}
}
close INFO;
if($name_old){
	$geneset{$name_old} = [@data];
}

open (INFO, $input_file2) or die $!;
open (OUT, ">$output_file") or die $!;
print OUT "# Column headers\n";
print OUT "# 1. Induced transcription factor (TF) gene symbol\n";
print OUT "# 2. Direction of gene expression change after TF induction\n";
print OUT "# 3. Target gene number ordered by increasing EPFP\n";
print OUT "# 4. Target gene symbol\n";
print OUT "# 5. Log10-ratio of target gene expression change after TF induction\n";
print OUT "# 6. All regulated targets (including surrogate ChIP-seq): EPFP\n";
print OUT "# 7. All regulated targets: Number of supporting ChIP-seq data sets\n";
print OUT "# 8. All regulated targets: genome position (center)\n";
print OUT "# 9. Direct regulated targets (no surrogate ChIP-seq): EPFP\n";
print OUT "# 10. Direct regulated targets: Number of supporting ChIP-seq data sets\n";
print OUT "# 11. Direct regulated targets: genome position (center)\n";
print OUT "# Abbreviations: N/S = no binding site or non-significant (EPFP>0.3); N/A = no ChIP-seq data\n";
print OUT "TF_name\tDirection\tNum\tSymbol\tLog10-ratio\tDirect-EPFP\tDirect-nData\tDirect-position\tAll-EPFP\tAll-nData\tAll-position\n";
my $output_table;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	if($line=~/^#/){ next; }
	my($name,$junk,@genes)=split(/\t/, $line);
	my ($TF,$dir) = split(/_/,$name);
	my $ref = $geneset{$name};
	my @data;
	my %hash;
	if($ref){
		@data = @$ref;
		my $N1 = @{$data[0]};
		for(my $i=0; $i<$N1; $i++){
			$hash{$data[0]->[$i]}=$i+1;
		}
	}
	my @data1;
	for(my $i=0; $i<4; $i++){
		$line = <INFO>;
		$line =~ s/\n$//;
		my($name1,$attr,@x)=split(/\t/, $line);
		push(@data1,\@x);
	}
	for(my $i=0; $i<@genes; $i++){
		my $num = $i+1;
		print OUT "$TF\t$dir\t$num\t$genes[$i]\t$data1[2]->[$i]\t$data1[0]->[$i]\t$data1[1]->[$i]\t$data1[3]->[$i]";
		if(!$direct{$TF}){
			print OUT "\tN/A\t0\tN/A\n";
			next;
		}
		my $i1 = $hash{$genes[$i]};
		if(!$i1){
			print OUT "\tN/S\t$direct{$TF}\tN/S\n";
			next;
		}
		$i1--;
		print OUT "\t$data[1]->[$i1]\t$data[2]->[$i1]\t$data[4]->[$i1]\n";
	}
}
close INFO;
close OUT;
exit(0);


