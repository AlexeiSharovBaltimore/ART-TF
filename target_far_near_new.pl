#!/usr/local/bin/perl

#target_far_near_new.pl TF_regul_targets_far_combined.txt TF_regul_targets_near_combined.txt regul_targets_table.txt
#target_far_near_new.pl TF_regul_targets_far_combined.txt TF_regul_targets_near_combined.txt regul_all_targets_genesets.txt -gall -anova anova-TFs_Nov2017_all.txt
#target_far_near_new.pl TF_regul_targets_add_far_combined.txt TF_regul_targets_add_near_combined.txt regul_all_targets_genesets_add.txt -gall -anova anova-TFs_Sep2018_all.txt -coord TF_targets_human_add_far2018.txt TF_targets_human_add_near2018.txt -indir indirect_list.txt

use strict;
my $usage = "$0 input_coord TSScoord output\n";
my $arg=0;
my $input_far = $ARGV[$arg++] or die $!;
my $input_near = $ARGV[$arg++] or die $!;
my $output_file = $ARGV[$arg++] or die $!;
my $option_geneset=0;
my $anova_file;
my $targets_file_far;
my $targets_file_near;
my $indirect_file;
my $output_file_response;
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-GDIR"){ $option_geneset=1; }
	elsif(uc($option) eq "-GALL"){ $option_geneset=2; }
	elsif(uc($option) eq "-GIND"){ $option_geneset=3; }
	elsif(uc($option) eq "-ANOVA"){ $anova_file = $ARGV[$arg++] or die $!; }
	elsif(uc($option) eq "-RESP"){ $output_file_response = $ARGV[$arg++] or die $!; }
	elsif(uc($option) eq "-INDIR"){ $indirect_file = $ARGV[$arg++] or die $!; }
	elsif(uc($option) eq "-COORD"){
		$targets_file_far = $ARGV[$arg++] or die $!;
		$targets_file_near = $ARGV[$arg++] or die $!;
	}
	else{ die "ERROR: Wrong option $option\n"; }
}

my @direction = ("up","down");
my $MISSING=-9999;
my %hashOriginal=("AHSA2P","AHSA2","ATP5F1A","ATP5A1","ATP5F1AP1","ATP5A1P1","ATP5F1AP2","ATP5A1P2","ATP5F1AP3","ATP5A1P3","ATP5F1AP7","ATP5A1P7","ATP5F1AP8","ATP5A1P8","ATP5F1B","ATP5B","ATP5F1C","ATP5C1","ATP5F1CP1","ATP5C1P1","ATP5F1D","ATP5D","ATP5F1E","ATP5E","ATP5F1EP2","ATP5EP2","ATP5PB","ATP5F1","ATP5PBP2","ATP5F1P2","ATP5PBP3","ATP5F1P3","ATP5PBP5","ATP5F1P5","ATP5PBP6","ATP5F1P6","ATP5PBP7","ATP5F1P7","ATP5MC1","ATP5G1","ATP5MC1P4","ATP5G1P4","ATP5MC1P5","ATP5G1P5","ATP5MC1P7","ATP5G1P7","ATP5MC2","ATP5G2","ATP5MC2P1","ATP5G2P1","ATP5MC2P3","ATP5G2P3","ATP5MC2P4","ATP5G2P4","ATP5MC3","ATP5G3","ATP5PD","ATP5H","ATP5PDP2","ATP5HP2","ATP5PDP3","ATP5HP3","ATP5PDP4","ATP5HP4","ATP5ME","ATP5I","ATP5PF","ATP5J","ATP5MF","ATP5J2","ATP5MFP4","ATP5J2P4","ATP5MFP5","ATP5J2P5","ATP5MF-PTCD1","ATP5J2-PTCD1","ATP5PFP1","ATP5JP1","ATP5MG","ATP5L","ATP5MGL","ATP5L2","ATP5MGP2","ATP5LP2","ATP5MGP3","ATP5LP3","ATP5MGP4","ATP5LP4","ATP5MGP5","ATP5LP5","ATP5PO","ATP5O","ATP5IF1","ATPIF1","DEPP1","C10orf10","ATP5MPL","C14orf2","CFAP97D1","C17orf105","RMC1","C18orf8","MIR29B2CHG","C1orf132","RHEX","C1orf186","ODR4","C1orf27","MIR1-1HG-AS1","C20orf166-AS1","SNORC","C2orf82","CAND1","CAND1.11","EFL1P1","EFTUD1P1","EP400P1","EP400NL","GUCY1A1","GUCY1A3","GUCY1B1","GUCY1B3","TRMT9B","KIAA1456","IQCN","KIAA1683","SMIM37","LINC00116","MIR570HG","LINC00969","ROCR","LINC02095","MIR4453HG","LINC02486","RNF227","LINC02581","SLC49A3","MFSD7","NUDT4B","NUDT4P1","PRSS46P","PRSS46","COP1","RFWD2","MTREX","SKIV2L2","TBXT","T","STIMATE","TMEM110","CLTRN","TMEM27","RXYLT1","TMEM5","MACO1","TMEM57","ATP5MD","USMG5","ATP5MDP1","USMG5P1","RAB6D","WTH3DI","ZRSR2P1","ZRSR1");

my @gene_symbols;
my %hashGene;
my %hashTargets;
my %hashCoord;
my %TF_table;
my %hashIndir;
if($indirect_file){
	open (INFO, $indirect_file) or die "file $indirect_file NOT found\n";
	while(my $line = <INFO>){
		chop $line;
		my ($TF0,@indirect)=split(/\t/, $line);
		foreach my $TF1 (@indirect){
			if(!$hashIndir{$TF1}){ $hashIndir{$TF1} = []; }
			push(@{$hashIndir{$TF1}},$TF0);
		}
	}
	close INFO;
}
if($targets_file_far){
	open (INFO, $targets_file_far) or die "file $targets_file_far NOT found\n";
	while(my $line = <INFO>){
		chop $line;
		if(!$line || $line =~ /^!|^</){ next; }
		my ($name,$descr,@symbols)=split(/\t/, $line);
		if(!$name){ next; }
		my ($TF,$junk) = split(/_/, $name);
		my $line1 = <INFO>;
		if($line1 !~ /^\tscore\t/){ print "ERR: No score\n"; exit(0); }
		chop $line1;
		my ($junk1,$junk2,@score)=split(/\t/, $line1);
		$line1 = <INFO>;
		if($line1 !~ /^\tposition\t/){ print "ERR: No position\n"; exit(0); }
		chop $line1;
		my ($junk3,$junk4,@coord)=split(/\t/, $line1);
		for(my $i=0; $i<@symbols; $i++){
			my $symbol = $symbols[$i];
			my $i1 = $hashGene{$symbol};
			if(!$i1){
				push(@gene_symbols,$symbol);
				$hashGene{$symbol} = @gene_symbols;
				$i1 = @gene_symbols;
			}
			my $ref = $hashCoord{"$TF,$i1"};
			if(!$ref){
				$hashCoord{"$TF,$i1"} = [$score[$i],$coord[$i]];
			}elsif($ref->[0] < $score[$i]){
				$ref->[0] = $score[$i];
				$ref->[1] = $coord[$i];
			}
		}
	}
	close INFO;
}
open (INFO, $input_far) or die "$input_far NOT found\n";
while(my $line = <INFO>){
	chop $line;
	my ($name,$Ngenes,@symbols)=split(/\t/, $line);
	if(!$name){ next; }
	my ($TF,$junk) = split(/_/, $name);
	my $line1 = <INFO>;
	chop $line1;
	my ($junk1,$junk2,@epfp)=split(/\t/, $line1);
	my $line2 = <INFO>;
	chop $line2;
	my ($junk3,$junk4,@nData)=split(/\t/, $line2);
	my @list;
	#my %hashNred;
	for(my $i=0; $i<@symbols; $i++){
		my $symbol = $symbols[$i];
		if($hashOriginal{$symbol}){ $symbol=$hashOriginal{$symbol}; }
		if($symbol eq $TF){ next; }
		#if($hashNred{$symbol}){ next; }
		#$hashNred{$symbol} = 1;
		my $i1 = $hashGene{$symbol};
		if(!$i1){
			push(@gene_symbols,$symbol);
			$hashGene{$symbol} = @gene_symbols;
			$i1 = @gene_symbols;
		}
		$hashTargets{$TF}->[$i1]=1000.001;
		my @TFset = ($TF);
		if($name =~ /_indir/ && $hashIndir{$TF}){ @TFset = @{$hashIndir{$TF}}; }
		my ($score,$coord) = (0,"");
		foreach my $TF1 (@TFset){
			my $ref1 = $hashCoord{"$TF1,$i1"};
			if($ref1 && $score < $ref1->[0]){
				($score,$coord) = @$ref1;
			}
		}
		if($targets_file_far && !$score){
			print "ERR: No coord far $TF $i1 $symbol\n";
		}
		push(@list,[$i1,$epfp[$i],$nData[$i],$coord]);
	}
	if(!@list){ next; }
	@list = sort {$a->[0]<=>$b->[0]} @list;
	#my $num = @list;
	#print "$name\n$num";
	my ($TF,$indir,$dir) = split(/_/,$name);
	if($name =~ /_up/){ $dir=0; }
	else{ $dir=1; }
	if($name =~ /_indir/){ $indir=1; }
	else{ $indir=0; }
	my $pos = $dir*4+$indir;
	$TF_table{$TF}->[$pos] = \@list;
	my $nn2 = @list;
}	
close INFO;
%hashCoord=();
if($targets_file_near){
	open (INFO, $targets_file_near) or die "file $targets_file_near NOT found\n";
	while(my $line = <INFO>){
		chop $line;
		if(!$line || $line =~ /^!|^</){ next; }
		my ($name,$descr,@symbols)=split(/\t/, $line);
		if(!$name){ next; }
		my ($TF,$junk) = split(/_/, $name);
		my $line1 = <INFO>;
		if($line1 !~ /^\tscore\t/){ print "ERR: No score\n"; exit(0); }
		chop $line1;
		my ($junk1,$junk2,@score)=split(/\t/, $line1);
		$line1 = <INFO>;
		if($line1 !~ /^\tposition\t/){ print "ERR: No position\n"; exit(0); }
		chop $line1;
		my ($junk,$junk,@coord)=split(/\t/, $line1);
		for(my $i=0; $i<@symbols; $i++){
			my $symbol = $symbols[$i];
			my $i1 = $hashGene{$symbol};
			if(!$i1){
				push(@gene_symbols,$symbol);
				$hashGene{$symbol} = @gene_symbols;
				$i1 = @gene_symbols;
			}
			my $ref = $hashCoord{"$TF,$i1"};
			if(!$ref){
				$hashCoord{"$TF,$i1"} = [$score[$i],$coord[$i]];
			}elsif($ref->[0] < $score[$i]){
				$ref->[0] = $score[$i];
				$ref->[1] = $coord[$i];
			}
		}
	}	
	close INFO;
}
open (INFO, $input_near) or die "$input_far NOT found\n";
while(my $line = <INFO>){
	chop $line;
	my ($name,$Ngenes,@symbols)=split(/\t/, $line);
	if(!$name){ next; }
	my ($TF,$junk) = split(/_/, $name);
	my $line1 = <INFO>;
	chop $line1;
	my ($junk1,$junk2,@epfp)=split(/\t/, $line1);
	$line1 = <INFO>;
	chop $line1;
	my ($junk3,$junk4,@nData)=split(/\t/, $line1);
	my @list;
	#my %hashNred;
	for(my $i=0; $i<@symbols; $i++){
		my $symbol = $symbols[$i];
		if($hashOriginal{$symbol}){ $symbol=$hashOriginal{$symbol}; }
		if($symbol eq $TF){ next; }
		#if($hashNred{$symbol}){ next; }
		#$hashNred{$symbol} = 1;
		my $i1 = $hashGene{$symbol};
		if(!$i1){
			push(@gene_symbols,$symbol);
			$hashGene{$symbol} = @gene_symbols;
			$i1 = @gene_symbols;
		}
		$hashTargets{$TF}->[$i1]=1000.001;
		my @TFset = ($TF);
		if($name =~ /_indir/ && $hashIndir{$TF}){ @TFset = @{$hashIndir{$TF}}; }
		my ($score,$coord) = (0,"");
		foreach my $TF1 (@TFset){
			my $ref1 = $hashCoord{"$TF1,$i1"};
			if($ref1 && $score < $ref1->[0]){
				($score,$coord) = @$ref1;
			}
		}
		if($targets_file_far && !$score){
			print "ERR: No coord near $TF $i1 $symbol\n";
		}
		push(@list,[$i1,$epfp[$i],$nData[$i],$coord]);
	}
	if(!@list){ next; }
	#my $num = @list;
	#print "$name\n$num";
	@list = sort {$a->[0]<=>$b->[0]} @list;
	my ($TF,$dir,$indir) = split(/_/,$name);
	if($name =~ /_up/){ $dir=0; }
	else{ $dir=1; }
	if($name =~ /_indir/){ $indir=1; }
	else{ $indir=0; }
	my $pos = 2+$dir*4+$indir;
	$TF_table{$TF}->[$pos] = \@list;
}	
close INFO;

if($option_geneset){
	print_geneset_file();
	exit(0);
}
open (OUT, ">$output_file") or die $!;
print OUT join("\t","TF","Far Up","Far Up all","Near Up","Near Up all","Far Down","Far Down all","Near Down","Near Down all","Overlap Up","Overlap Up all","Overlap Down","Overlap Down all")."\n";
foreach my $TF (sort keys %TF_table){
	print "$TF\n";
	my $ref = $TF_table{$TF};
	my @Ngenes;
	for(my $dir=0; $dir<2; $dir++){
		for(my $indir=0; $indir<2; $indir++){
			my ($pos1,$pos2,$pos3)=($dir*4+$indir,2+$dir*4+$indir,8+$dir*2+$indir);
			my $ref1 = $ref->[$pos1];
			my $ref2 = $ref->[$pos2];
			if($indir && $ref->[$pos1-1]){
				if(!$ref1){ $ref1 = $ref->[$pos1-1]; }
				else{
					push(@$ref1,@{$ref->[$pos1-1]});
					@$ref1 = sort {$a->[0]<=>$b->[0]} @$ref1;
				}
			}
			if($indir && $ref->[$pos2-1]){
				if(!$ref2){ $ref2 = $ref->[$pos2-1]; }
				else{
					push(@$ref2,@{$ref->[$pos2-1]});
					@$ref2 = sort {$a->[0]<=>$b->[0]} @$ref2;
				}
			}
			if (!$ref1){ $Ngenes[$pos1]=0; }
			else{ $Ngenes[$pos1]=@$ref1; }
			if (!$ref2){ $Ngenes[$pos2]=0; }
			else{ $Ngenes[$pos2]=@$ref2; }
			$Ngenes[$pos3]=0;
			if($ref1 && $ref2){
				my ($n12,$k1,$k2)=(0,0,0);
				while($k1<$Ngenes[$pos1] && $k2<$Ngenes[$pos2]){
					if($ref1->[$k1]->[0]==$ref2->[$k2]->[0]){ $k1++; $k2++; $n12++; }
					elsif($ref1->[$k1]->[0] < $ref2->[$k2]->[0]){ $k1++; }
					elsif($ref1->[$k1]->[0] > $ref2->[$k2]->[0]){ $k2++; }
				}
				$Ngenes[$pos3] = $n12;
			}
		}
	}
	print OUT "$TF\t".join("\t",@Ngenes)."\n";
}
close OUT;
exit(0);

#**********************************
sub  print_geneset_file
#**********************************
{
if($anova_file){
	open (INFO, $anova_file) or die $!;
	my $line = <INFO>;
	chop $line;
	my ($junk1,$junk2,$junk3,$junk4,@headers);
	my %hashTFs;
	my $ii=0;
	($junk1,$junk2,@headers) = split(/\t/,$line);
	while($ii<@headers && $headers[$ii] !~ /^Var/){ $ii++; }
	for(my $i=0; $i<$ii; $i++){
		my $TF1 = $headers[$i];
		my ($TF,$junk6) = split(/_| \(/,$TF1);
		$headers[$i] = $TF;
	}
	while(my $line = <INFO>){
		chop $line;
		my ($symbol,$junk3,@data) = split(/\t/, $line);
		if($hashOriginal{$symbol}){ $symbol=$hashOriginal{$symbol}; }
		my $i1 = $hashGene{$symbol};
		if(!$i1){ next; }
		splice(@data,$ii);
		my $median = median(\@data);
		for(my $i=0; $i<$ii; $i++){
			if($data[$i]==$MISSING){ next; }
			my $ref = $hashTargets{$headers[$i]};
			if(!$ref){ next; }
			if($ref->[$i1]>1000){
				$ref->[$i1] = int(10000*($data[$i]-$median))/10000;
			}
		}
	}
	close INFO;
}
if($output_file_response){ open (RESP, ">$output_file_response"); }
open (OUT, ">$output_file") or die $!;
foreach my $TF (sort keys %TF_table){
	my $ref = $TF_table{$TF};
	my @Ngenes;
	for(my $dir=0; $dir<2; $dir++){
		my ($pos1,$pos2)=($dir*4,2+$dir*4);
		my $ref1 = $ref->[$pos1];
		my $ref2 = $ref->[$pos2];
		if(!$ref1){ $ref->[$pos1] = []; $ref1 = $ref->[$pos1]; }
		if($option_geneset==2 && $ref->[$pos1+1]){ push(@$ref1,@{$ref->[$pos1+1]}); }
		if($ref2){ push(@$ref1,@$ref2); }
		if($option_geneset==2 && $ref->[$pos2+1]){ push(@$ref1,@{$ref->[$pos2+1]}); }
		if($option_geneset==3){
			@$ref1 = ();
			if($ref->[$pos1+1]){ push(@$ref1, @{$ref->[$pos1+1]}); }
			if($ref->[$pos2+1]){ push(@$ref1, @{$ref->[$pos2+1]}); }
		}
		@$ref1 = sort {$a->[0]<=>$b->[0]} @$ref1;
		for(my $i=0; $i<@$ref1; $i++){
			my ($s1,$epfp1,$nData1,$coord1) = @{$ref1->[$i]};
			for(my $j=$i+1; $j<@$ref1; $j++){
				my ($s2,$epfp2,$nData2,$coord2) = @{$ref1->[$j]};
				if($s2 == $s1){
					if($epfp1 > $epfp2){
						$ref1->[$i]->[1]=$epfp2;
					}
					if($coord1 && $coord2){
						$ref1->[$i]->[3] .= ",$coord2";
					}elsif($coord2){
						$ref1->[$i]->[3] .= $coord2;
					}
					$ref1->[$i]->[2] += $nData2;
					splice(@$ref1,$j--,1);
				}
			}
		}
		my @string1; my @string2; my @string3; my @string4; my @string5;
		my $target_resp=0;
		foreach my $ref4 (sort {$a->[1]<=>$b->[1]} @$ref1){
			my ($is,$epfp,$nData,$coord) = @$ref4; 
			$epfp = int(100000*$epfp)/100000;
			push(@string1,$gene_symbols[$is-1]);
			push(@string2,$epfp);
			push(@string3,$nData);
			if($anova_file){
				my $x3 = $hashTargets{$TF}->[$is];
				if($x3>1000){
					print "$gene_symbols[$is-1]\n";
					$x3 = 0.301;
					if($dir){ $x3 = -0.301; }
				}
				push(@string4,$x3);
				$target_resp++;
			}
			if($targets_file_far){
				push(@string5,$coord);
			}
		}
		my $nn1 = @string1;
		if(!$nn1){ next; }
		my $dirWord = "_$direction[$dir]";
		if($output_file_response){
			print RESP "$TF"."_$direction[$dir]\t$nn1\t$target_resp\n";
		}
		print OUT "$TF$dirWord\t$nn1\t".join("\t",@string1)."\n";
		print OUT "\tEPFP\t".join("\t",@string2)."\n";
		print OUT "\tnData\t".join("\t",@string3)."\n";
		if($anova_file){
			print OUT "\tlogratio\t".join("\t",@string4)."\n";
		}
		if($targets_file_far){
			print OUT "\tposition\t".join("\t",@string5)."\n";
		}
	}
}
close OUT;
if($output_file_response){ close RESP; }
return;
}

#**************************
sub   median
#**************************
{
my $xr = shift;
if(!$xr){ return 0; }
my $n = @$xr;
if(!$n){ return 0; }
my $median = 0;
my @sorted = sort {$a<=>$b} @$xr;
while($sorted[0]==$MISSING){
	shift(@sorted);
	$n--;
}
my $i = $n/2;
if($i > int($i)){
	$i = int($i);
	$median = $sorted[$i];
}else{
	$median = ($sorted[$i] + $sorted[$i-1])/2;
}
return $median;
}

