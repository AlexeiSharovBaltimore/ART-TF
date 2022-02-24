#!/usr/local/bin/perl

#output2targets.pl CREST-TF_regulated_targets_near_EPFP03_new TF_regul_targets_near.txt -ind indirect_list.txt

use strict;
my $usage = "$0 input_coord TSScoord output\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $!;
my $output_file = $ARGV[$arg++] or die $!;
my $indirect_file;
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-IND"){ $indirect_file = $ARGV[$arg++] or die $usage; }
	else{ die "ERROR: Wrong option $option\n"; }
}
my @direction = ("up","down");
my @removeTargets=("ARNT2_01","ATF2_52","ATF3_14","DLX4_01","DLX6_01","DMRT2_01","DUX4_50","GLI1_01","HHEX_01","HNF4A_05","IRX2_01","KLF15_01","KLF5_07","KLF5_08","NKX2-5_01","PATZ1_01","SALL4_03","SOX9_02","TFE3_01","TP53_50","CHD2_01","EZH2_01","EZH2_02","EZH2_03","H3K27ac_01","H3K27me3_01","H3K9ac_01","H3K9me3_01","H4K20me1_01","HAND1_01","HAND1_02","HAND1_03","HAND2_01","HDAC1_01","HOXB13_01","MBD2_01","MBD3L2_01","MBD3L2_02","MBD3L2_03","MED1_01","MED12_01","MSX1_01","NR5A2_01","NR5A2_02","P300_01","P300_02","P300_03","POL2_01","POL2_02","POL2_03","POL2_04","POL2_05","POL2_06","POL2_07","POL3_01","POL3_02","SP1_01","SP1_02","SP1_03","SP1_04","SP1_05","STAT3_02","STAT3_04","STAT3_05","STAT3_07","STAT3_08","STAT3_10","STAT3_11","STAT3_12","STAT1_04","TBP_01");

my %hashIndirect;
if($indirect_file){
	open (INFO, $indirect_file) or die "$indirect_file NOT found\n";
	while(my $line = <INFO>){
		if($line =~ /^#/){ next; }
		chop $line;
		my ($TF2,@TF1list) = split(/\t/, $line);
		for(my $i=0; $i<@TF1list; $i++){
			my $TF1 = $TF1list[$i];
			$hashIndirect{$TF2}->{$TF1}=1;
		}
	}
	close INFO;
}
my %removeTargets;
foreach my $target (@removeTargets){
	$removeTargets{$target}=1;
}

my @targetList;
my %hashTarget;
my @DATA1;
my @DATAindirect;
my @geneset;
open (INFO, $input_file) or die "$input_file NOT found\n";
while(my $line = <INFO>){ if($line =~ /^!Matrix4_start/){ last; } }
my $line = <INFO>;
chop $line;
my ($junk,@TFname)=split(/\t/, $line);
my %hash;
while(my $line = <INFO>){
	chop $line;
	if($line =~ /^!/){ last; }
	my ($targets,@data)=split(/\t/, $line);
	my @items = split(/_/,$targets);
	if(@items==3){
		if($targets =~ /^RBBP5_KDM5/){ $targets=$items[0]."_$items[2]"; }
		else{ $targets=$items[1]."_$items[2]"; }
	}
	if($removeTargets{$targets}){ next; }
	if($hash{$targets}){ print "ERR: Duplicated $targets\n"; exit(0); }
	$hash{$targets}=1;
	my $itarget = $hashTarget{$targets};
	if(!$itarget){
		push(@targetList,$targets);
		$hashTarget{$targets} = @targetList;
		$itarget = @targetList;
	}
	my $target1 = $targets;
	$target1 =~ s/_.+$//;
	my $ref = $hashIndirect{$target1};
	for(my $i=0; $i<@data; $i++){
		if(!$data[$i]){ next; }
		my $name = $TFname[$i];
		if($TFname[$i] eq $target1){
			$DATA1[$itarget-1]->[0] = $data[$i];
		}elsif($ref && $ref->{$name}){
			$DATAindirect[$itarget-1]->{$name}->[0] = $data[$i];
		}
	}
}

while(my $line = <INFO>){ if($line =~ /^!Matrix5_start/){ last; } }
$line = <INFO>;
chop $line;
my %hash=();
while(my $line = <INFO>){
	chop $line;
	if($line =~ /^!/){ last; }
	my ($targets,@data)=split(/\t/, $line);
	my @items = split(/_/,$targets);
	if(@items==3){ $targets=$items[1]."_$items[2]"; }
	if($removeTargets{$targets}){ next; }
	if($hash{$targets}){ print "ERR: Duplicated $targets\n"; exit(0); }
	$hash{$targets}=1;
	my $itarget = $hashTarget{$targets};
	my $target1 = $targets;
	$target1 =~ s/_.+$//;
	my $ref = $hashIndirect{$target1};
	for(my $i=0; $i<@data; $i++){
		if(!$data[$i]){ next; }
		my $name = $TFname[$i];
		if($TFname[$i] eq $target1){
			$DATA1[$itarget-1]->[1] = $data[$i];
		}elsif($ref && $ref->{$name}){
			$DATAindirect[$itarget-1]->{$name}->[1] = $data[$i];
		}
	}
}
while(my $line = <INFO>){ if($line =~ /^!Database_start/){ last; } }
while(my $line = <INFO>){
	chop $line;
	if($line =~ /^>/){
		my $ID = $line;
		$ID =~ s/^>//;
		$line = <INFO>;
		chop $line;
		$line =~ s/^\S+\t//;
		my @genes=split(/,/,$line);
		$geneset[$ID]->[0] = \@genes;
		$line = <INFO>;
		chop $line;
		$line =~ s/^\S+\t//;
		my @epfp=split(/,/,$line);
		$geneset[$ID]->[1] = \@epfp;
	}
}
close INFO;
open (OUT, ">$output_file") or die $!;
for(my $it=0; $it<@targetList; $it++){
	my $ref = $DATA1[$it];
	if($ref){
		my $TF = $targetList[$it];
		$TF =~ s/_.+$//;
		for(my $idir=0; $idir<2; $idir++){
			my $ID = $ref->[$idir];
			if(!$ID){ next; }
			print OUT "$targetList[$it]//$TF"."_$direction[$idir]\t\t"; 
			print OUT join("\t",@{$geneset[$ID]->[0]})."\n";
			print OUT "\tEPFP\t"; 
			print OUT join("\t",@{$geneset[$ID]->[1]})."\n";
		}
	}
	$ref = $DATAindirect[$it];
	if(!$ref){ next; }
	foreach my $TF (keys %$ref){
		my $ref1 = $ref->{$TF};
		for(my $idir=0; $idir<2; $idir++){
			my $ID = $ref1->[$idir];
			if(!$ID){ next; }
			print OUT "$targetList[$it]//$TF"."_$direction[$idir]\t\t"; 
			print OUT join("\t",@{$geneset[$ID]->[0]})."\n";
			print OUT "\tEPFP\t"; 
			print OUT join("\t",@{$geneset[$ID]->[1]})."\n";
		}

	}
}
close OUT;
exit(0);


