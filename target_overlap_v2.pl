#!/usr/local/bin/perl

#target_overlap.pl TF_regul_targets_add_far2021.txt ndata_chip.txt TF_regul_targets_add_far_combined2022.txt -g

use strict;
my $usage = "$0 input_coord TSScoord output\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $!;
my $ndata_file = $ARGV[$arg++] or die $!;
my $output_file = $ARGV[$arg++] or die $!;
my $option_geneset=0;
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-G"){ $option_geneset=1; }
	else{ die "ERROR: Wrong option $option\n"; }
}
my @direction = ("up","down");

my %ndata;
open (INFO, $ndata_file) or die "$input_file NOT found\n";
while(my $line = <INFO>){
	chop $line;
	my ($TF,@nnn)=split(/\t/, $line);
	$ndata{$TF} = \@nnn;
}
close INFO;
my @gene_symbols;
my %hashGene;
my %hashGeneset;
my %hashTF;
my %hashTarget;
my %group2compare;
open (INFO, $input_file) or die "$input_file NOT found\n";
while(my $line = <INFO>){
	chop $line;
	my ($name,$descr,@symbols)=split(/\t/, $line);
	if(!$name){ next; }
	my ($TF,$junk) = split(/_/, $name);
	my $line1 = <INFO>;
	chop $line1;
	my ($junk1,$junk2,@epfp)=split(/\t/, $line1);
	my @list;
	for(my $i=0; $i<@symbols; $i++){
		my $symbol = $symbols[$i];
		if($symbol eq $TF){ next; }
		my $i1 = $hashGene{$symbol};
		if(!$i1){
			push(@gene_symbols,$symbol);
			$hashGene{$symbol} = @gene_symbols;
			$i1 = @gene_symbols;
		}
		push(@list,[$i1,$epfp[$i]]);
	}
	if(!@list){ next; }
	@list = sort {$a->[0]<=>$b->[0]} @list;
	#print "$name\n@list";

	$hashGeneset{$name} = \@list;
	my ($target,$ending) = split(/\/\//,$name);
	my ($TF,$dir) = split(/_/,$ending);
	if($dir eq "up"){ $dir=0; }
	else{ $dir=1; }
	push(@{$group2compare{$TF}->[$dir]}, $target);
}	
close INFO;

open (OUT, ">$output_file") or die $!;
foreach my $TF (sort keys %group2compare){
	for(my $dir=0; $dir<2; $dir++){
		my $ref = $group2compare{$TF}->[$dir];
		if (!$ref){ next; }
		my @targets = sort @$ref;
		my $template = $TF."_\\d+";
		my $i1=0;
		my @native;
		while($i1<@targets && $targets[$i1] !~ /^$template/){ $i1++; }
		if($i1<@targets){
			my $i2=$i1;
			while($i2<@targets && $targets[$i2] =~ /^$template/){ $i2++; }
			@native = splice(@targets,$i1,$i2-$i1);
			unshift(@targets,@native);
		}
		my %hashIndirect;
		my %hashDirect;
		for(my $i=0; $i<@targets; $i++){
			my $name1 = $targets[$i]."//".$TF."_".$direction[$dir];
			my $ref1 = $hashGeneset{$name1};
			if(!$ref1){ next; }
			my @list1 = @$ref1;
			my $N1 = @list1;
			foreach my $ref2 (@list1){
				my($is,$epfp)=@$ref2;
				my $ref3;
				if($i < @native){ $ref3 = $hashDirect{$is}; }
				else{ $ref3 = $hashIndirect{$is}; }
				if($ref3){
					$ref3->[0]++;
					if($ref3->[1]>$epfp){ $ref3->[1] = $epfp; }
				}else{
					if($i < @native){ $hashDirect{$is} = [1,$epfp]; }
					else{ $hashIndirect{$is} = [1,$epfp]; }
				}
			}
		}
		my @listDirect;
		foreach my $key (keys %hashDirect){
			my ($x1,$x2) = @{$hashDirect{$key}};
			if($ndata{$TF}->[1] <3 || $x1>1){
				push(@listDirect,[$key,$x2,$x1]);
				$hashDirect{$key}->[0] = 100000;
			}
		}
		if(@listDirect && $option_geneset){
			my @string1; my @string2; my @string3;
			foreach my $ref2 (sort {$a->[1]<=>$b->[1]} @listDirect){
				my ($is,$x2,$x1) = @$ref2; $x2 = int(100000*$x2)/100000;
				push(@string1,$gene_symbols[$is-1]);
				push(@string2,$x2);
				push(@string3,$x1);
			}
			my $nn1 = @listDirect;
			print OUT "$TF"."_$direction[$dir]\t$nn1\t".join("\t",@string1)."\n";
			print OUT "\tEPFP\t".join("\t",@string2)."\n";
			print OUT "\tnData\t".join("\t",@string3)."\n";
		}
		my @listIndirect;
		foreach my $key (keys %hashIndirect){
			my ($x1,$x2) = @{$hashIndirect{$key}};
			my $ref3 = $hashDirect{$key};
			if($ref3){
				$x1 += $ref3->[0]; $x2 += $ref3->[1];
			}
			if($ndata{$TF}->[0] <3 || $x1>1 && $x1<100000){
				push(@listIndirect,[$key,$x2,$x1]);
			}
		}
		if(@listIndirect && $option_geneset){
			my @string1; my @string2; my @string3;
			foreach my $ref2 (sort {$a->[1]<=>$b->[1]} @listIndirect){
				my ($is,$x2,$x1) = @$ref2; $x2 = int(100000*$x2)/100000;
				push(@string1,$gene_symbols[$is-1]);
				push(@string2,$x2);
				push(@string3,$x1);
			}
			my $nn1 = @listIndirect;
			print OUT "$TF"."_indir_$direction[$dir]\t$nn1\t".join("\t",@string1)."\n";
			print OUT "\tEPFP\t".join("\t",@string2)."\n";
			print OUT "\tnData\t".join("\t",@string3)."\n";
		}
		if($option_geneset){ next; }

		my @M;
		for(my $i=0; $i<@targets; $i++){
			my $name1 = $targets[$i]."//".$TF."_".$direction[$dir];
			my $ref1 = $hashGeneset{$name1};
			if(!$ref1){ next; }
			my @list1 = @$ref1;
			my $N1 = @list1;
			$M[$i]->[$i] = $N1;
			for(my $j=$i+1; $j<@targets; $j++){
				my $name2 = $targets[$j]."//".$TF."_".$direction[$dir];
				my $ref2 = $hashGeneset{$name2};
				if(!$ref2){ next; }
				my @list2 = @$ref2;
				my $N2 = @list2;
				my ($n12,$k1,$k2)=(0,0,0);
				while($k1<$N1 && $k2<$N2){
					if($list1[$k1]->[0]==$list2[$k2]->[0]){ $k1++; $k2++; $n12++; }
					elsif($list1[$k1]->[0] < $list2[$k2]->[0]){ $k1++; }
					elsif($list1[$k1]->[0] > $list2[$k2]->[0]){ $k2++; }
				}
				$M[$i]->[$j] = $n12;
				$M[$j]->[$i] = $n12;
			}
		}
		if(!@M){ last; }
		print OUT "$TF"."_$direction[$dir]\n";
		print OUT "Targets\t".join("\t",@targets)."\n";
		for(my $i=0; $i<@targets; $i++){
			print OUT "$targets[$i]\t".join("\t",@{$M[$i]})."\n";
		}
		print OUT "\n";
	}
}
close OUT;
exit(0);


