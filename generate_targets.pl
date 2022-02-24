#!/usr/local/bin/perl

#generate_targets.pl peaks_combined_single.txt targets_human_TFs_far.txt -o 1
#generate_targets.pl peaks_combined_single.txt targets_human_TFs_near.txt -o 2
#generate_targets.pl peaks_combined_single.txt targets_human_TFs_near_upstrm.txt -o 3
#generate_targets.pl peaks_combined_single.txt targets_human_TFs_near_dwnstrm.txt -o 4

use strict;
my $usage = "$0 input_coord TSScoord output\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $usage;
my $output_file = $ARGV[$arg++] or die $usage;
my $score_thresh = 0;
my $option_TFBS = 0;
my $remove_TF = 0;
my $add_binding_sites = 0;
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-TF"){ $remove_TF = 1; }
	elsif(uc($option) eq "-O"){ $option_TFBS = $ARGV[$arg++] or die $usage; }
	elsif(uc($option) eq "-ADD"){ $add_binding_sites = 1; }
	else{ die "ERROR: Wrong option $option\n"; }
}

my $OPTION_ALL=0;
my $OPTION_FAR=1;
my $OPTION_NEAR=2;
my $OPTION_NEAR_UPSTREAM=3;
my $OPTION_NEAR_DOWNSTREAM=4;

my %hashTargets;
my ($count,$TFold)=(0,"");
open (INFO, $input_file) or die $!;
my $line = <INFO>;
my $count=0;
my @peaks;
while(my $line = <INFO>){
	chop $line;
	my ($TF,$chr,$pos,$score,$Nsupport,$binCode,$symbols,$distances)=split(/\t/, $line);
	if($chr !~ /^chr/){ last; }
	if(!$symbols){ next; }
	if($TFold && $TF ne $TFold){
		print "$TFold\n";
		process_peaks($TFold,\@peaks);
		@peaks=();
	}
	push(@peaks,[$chr,$pos,$score,$Nsupport,$binCode,$symbols,$distances]);
	$TFold=$TF;
}
if($TFold){
	print "$TFold\n";
	process_peaks($TFold,\@peaks);
}
close INFO;

my $count1=1;
open (OUT, ">$output_file") or die $!;
foreach my $TF (sort keys %hashTargets){
	my $ref = $hashTargets{$TF};
	my @sorted = sort {$ref->{$b}<=>$ref->{$a}} keys %$ref;
	my $maxScore = $ref->{$sorted[0]};
	if(@sorted > 1000){ $maxScore = $ref->{$sorted[50]}; }
	elsif(@sorted > 500){ $maxScore = $ref->{$sorted[25]}; }
	elsif(@sorted > 200){ $maxScore = $ref->{$sorted[10]}; }
	elsif(@sorted > 100){ $maxScore = $ref->{$sorted[5]}; }
	elsif(@sorted > 50){ $maxScore = $ref->{$sorted[3]}; }
	elsif(@sorted > 20){ $maxScore = $ref->{$sorted[2]}; }
	my $count=0;
	if(@sorted >= 10){
		print OUT "$TF\t$TF";
		my $text1 = "\tscore";
		for(my $i=0; $i<@sorted; $i++){
			my $symbol=$sorted[$i];
			my $score = int(100000*$ref->{$symbol}/$maxScore)/1000;
			if($score < 0.01 || $i>=5000){ last; }
			print OUT "\t$symbol";
			$text1 .= "\t$score";
			$count++;
		}
		print OUT "\n$text1\n";
	}
	$maxScore = int(1000*$maxScore)/1000;
	print "$TF\t$maxScore\t$count1\n";
	$count1++;
}
close OUT;
exit(0);

#**********************************
sub process_peaks
#**********************************
{
my $TF = shift;
my $peak_ref = shift;

my ($TF1,$junk) = split(/_/,$TF);
my $n = @$peak_ref;
my @sorted = sort {$peak_ref->[$b]->[2]<=>$peak_ref->[$a]->[2]} 0..($n-1);
my %targets;
for(my $ipeak=0; $ipeak<$n && $ipeak<30000; $ipeak++){
	my ($chr,$pos,$score,$Nsupport,$binCode,$symbols,$distances) = @{$peak_ref->[$sorted[$ipeak]]};
	my @symbol = split(/,/,$symbols);
	my @dist = split(/,/,$distances);
	my @near;
	for(my $i=0; $i<@symbol; $i++){
		if($dist[$i]>=-500 && $dist[$i]<=500){
			if($dist[$i]>=0){ $near[$i]=1; }
			else{ $near[$i]=2; }
		}
	}
	for(my $i=0; $i<@symbol; $i++){
		my $symbol1 = $symbol[$i];
		if($remove_TF && $symbol1 eq $TF1){ next; }
		if($option_TFBS==$OPTION_FAR && @near){ next; }
		if($option_TFBS > $OPTION_FAR && !$near[$i]){ next; }
		if($option_TFBS==$OPTION_NEAR_UPSTREAM && $near[$i]!=1){ next; }
		if($option_TFBS==$OPTION_NEAR_DOWNSTREAM && $near[$i]!=2){ next; }
		my $d = abs($dist[$i])/1000;
		if($d<1){ $d=1; }
		if(!$targets{$symbol1}){
			$targets{$symbol1} = $score/$d;
		}else{
			if($add_binding_sites){
				$targets{$symbol1} += $score/$d;
			}elsif($targets{$symbol1} < $score/$d){
				$targets{$symbol1} = $score/$d;
			}
		}
	}
}
$hashTargets{$TF} = \%targets;
return;
}

