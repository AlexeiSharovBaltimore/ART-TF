#!/usr/local/bin/perl

#parse_targets.pl TF_regulated_targets_EPFP03_add_far.txt TF_regulated_targets_EPFP03_add_near.txt target_scores_add.txt -ind indirect_list.txt

use strict;
my $usage = "$0 input_coord TSScoord output\n";
my $arg=0;
my $input_file1 = $ARGV[$arg++] or die $!;
my $input_file2 = $ARGV[$arg++] or die $!;
my $output_file = $ARGV[$arg++] or die $!;
my $indirect_file;
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-IND"){ $indirect_file = $ARGV[$arg++] or die $usage; }
	else{ die "ERROR: Wrong option $option\n"; }
}

my @location = ("far","near");
my @direction = ("up","dwn");
my @targetList;
my %hashTarget;
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

my @DATA1;
my @DATAindirect;
my @Ngenes;
for(my $iloc=0; $iloc<2; $iloc++){
	my $input_file = $input_file1;
	if($iloc==1){ $input_file = $input_file2; }
	print "$input_file\n";
	if(!open (INFO, $input_file)){
		print "NOT found, skipping\n";
		next;
	}
	while(my $line = <INFO>){ if($line =~ /^!Matrix1_start/){ last; } }
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
			my $name = $TFname[$i];
			if($TFname[$i] eq $target1){  #direct targets
				my $x = $data[$i]; if($x<0){ $x=0; }
				$DATA1[$itarget-1]->[$iloc] = $x;
			}elsif($ref && $ref->{$name}){  #indirect targets
				my $x = $data[$i]; if($x<0){ $x=0; }
				$DATAindirect[$itarget-1]->{$name}->[$iloc] = $x;
			}
		}
	}
	while(my $line = <INFO>){ if($line =~ /^!Matrix2_start/){ last; } }
	my $line = <INFO>;
	chop $line;
	my ($junk,@TFname)=split(/\t/, $line);
	%hash=();
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
		my $target1 = $targets;
		$target1 =~ s/_.+$//;
		my $ref = $hashIndirect{$target1};
		for(my $i=0; $i<@data; $i++){
			my $name = $TFname[$i];
			if($TFname[$i] eq $target1){
				my $x = $data[$i]; if($x<0){ $x=0; }
				$DATA1[$itarget-1]->[2+$iloc] = $x;
			}elsif($ref && $ref->{$name}){
				my $x = $data[$i]; if($x<0){ $x=0; }
				$DATAindirect[$itarget-1]->{$name}->[2+$iloc] = $x;
			}
		}
	}
	while(my $line = <INFO>){ if($line =~ /^!Matrix4_start/){ last; } }
	my $line = <INFO>;
	chop $line;
	my ($junk,@TFname)=split(/\t/, $line);
	%hash=();
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
		my $target1 = $targets;
		$target1 =~ s/_.+$//;
		my $ref = $hashIndirect{$target1};
		for(my $i=0; $i<@data; $i++){
			if(!$data[$i]){ next; }
			my $name = $TFname[$i];
			if($TFname[$i] eq $target1){
				$DATA1[$itarget-1]->[4+$iloc] = $data[$i];
			}elsif($ref && $ref->{$name}){
				$DATAindirect[$itarget-1]->{$name}->[4+$iloc] = $data[$i];
			}
		}
	}
	while(my $line = <INFO>){ if($line =~ /^!Matrix5_start/){ last; } }
	my $line = <INFO>;
	chop $line;
	my ($junk,@TFname)=split(/\t/, $line);
	%hash=();
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
		my $target1 = $targets;
		$target1 =~ s/_.+$//;
		my $ref = $hashIndirect{$target1};
		for(my $i=0; $i<@data; $i++){
			if(!$data[$i]){ next; }
			my $name = $TFname[$i];
			if($TFname[$i] eq $target1){
				$DATA1[$itarget-1]->[6+$iloc] = $data[$i];
			}elsif($ref && $ref->{$name}){
				$DATAindirect[$itarget-1]->{$name}->[6+$iloc] = $data[$i];
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
			my @genes=split(/,/,$line);
			$Ngenes[$iloc]->[$ID] = @genes;
			$line = <INFO>;
		}
	}
}
close INFO;
open (OUT, ">$output_file") or die $!;
print OUT "Target\tTF\tDirect\tZ-far-up\tZ-near-up\tZ-far-down\tZ-near-down\tNgenes-far-up\tNgenes-near-up\tNgenes-far-down\tNgenes-near-down\n";
for(my $it=0; $it<@targetList; $it++){
	my $TF = $targetList[$it];
	$TF =~ s/_.+$//;
	my $ref = $DATA1[$it];
	if($ref){
		my $text ="$targetList[$it]\t$TF\t1"; 
		my $found = 0;
		for(my $idir=0; $idir<2; $idir++){
			for(my $iloc=0; $iloc<2; $iloc++){
				my $pos = $idir*2+$iloc;
				my $x = $ref->[$pos];
				if($x){ $found=1; }
				$text .= "\t$x";
			}
		}
		for(my $idir=0; $idir<2; $idir++){
			for(my $iloc=0; $iloc<2; $iloc++){
				my $pos = 4+$idir*2+$iloc;
				my $x = $ref->[$pos];
				if($x){ $x=$Ngenes[$iloc]->[$x]; }
				$text .= "\t$x";
			}
		}
		if($found){ print OUT "$text\n"; }
	}
	$ref = $DATAindirect[$it];
	if(!$ref){ next; }
	foreach my $TF (keys %$ref){
		my $text ="$targetList[$it]\t$TF\t0"; 
		my $found = 0;
		my $ref1 = $ref->{$TF};
		for(my $idir=0; $idir<2; $idir++){
			for(my $iloc=0; $iloc<2; $iloc++){
				my $pos = $idir*2+$iloc;
				my $x = $ref1->[$pos];
				if($x){ $found=1; }
				$text .= "\t$x";
			}
		}
		for(my $idir=0; $idir<2; $idir++){
			for(my $iloc=0; $iloc<2; $iloc++){
				my $pos = 4+$idir*2+$iloc;
				my $x = $ref1->[$pos];
				if($x){ $x=$Ngenes[$iloc]->[$x]; }
				$text .= "\t$x";
			}
		}
		if($found){ print OUT "$text\n"; }
	}
}
close OUT;
exit(0);


