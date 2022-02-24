#!/usr/local/bin/perl

#combine_peaks_TFs.pl peaks/file_list.txt peaks_combined.txt -TSS TSS_refseq.txt
use strict;

my $usage = "$0 input_file\n";
my $arg = 0;
my $inputFile = $ARGV[$arg++] or die $usage;
my $outputFile = $ARGV[$arg++] or die $usage;
my $TSS_file = 0;
while(my $option = $ARGV[$arg++]){
	if(uc($option) eq "-TSS"){ $TSS_file = $ARGV[$arg++] or die $usage; }
	else{ die "ERROR: Wrong option $option\n"; }
}
my $delta = 500;
my $maxDist = 100000;
my $MISSING=-9999;

my %hashTSS;
if($TSS_file){
	open (INFO, $TSS_file) or die $!;
	my $line = <INFO>;
	while(my $line = <INFO>){
		$line =~ s/\n$//;
		my($symbol,$chr,$strand,$TSS)=split(/\t/, $line);
		push(@{$hashTSS{$chr}}, [$TSS,$strand,$symbol]);
	}
	close INFO;
	foreach my $chr (sort keys %hashTSS){
		@{$hashTSS{$chr}} = sort {$a->[0]<=>$b->[0]} @{$hashTSS{$chr}};
	}
}

my %files;
open (INFO, $inputFile) or die $!;
while(my $line = <INFO>){
	$line =~ s/\n$//;
	my ($TF,$file)=split(/\t/, $line);
	if($file){
		my ($id,@items) = split(/_/,$file);
		push(@{$files{$TF}},[$file,$id]);
	}
}
close(INFO);
print "TSS loaded\n";
my @report;
open (OUT, ">$outputFile") or die $!;
print OUT "TF\tchr\tpos\tscore\tNfiles\tbinCode\n";
foreach my $TF (sort keys %files){
	my $ref = $files{$TF};
	my $nfiles = @$ref;
	my %coord;
	my @peaks;
	my @scoreAvg;
	my @numPeaks;
	for(my $ifile=0; $ifile<$nfiles; $ifile++){
		my $found=0;
		my ($sum,$nn)=(0,0);
		my $filename = $ref->[$ifile]->[0];
		my $id = $ref->[$ifile]->[1];
		print "$TF $filename\n";
		if(!open (INFO, "$filename")){ print "ERR: $filename not found"; next; }
		while(my $line = <INFO>){
			$line =~ s/\n$//;
			my ($chr,$start,$end,$score)=split(/\t/, $line);
			if($end<$start){ die "Err: $end < $start in $filename\n"; }
			my $pos = int(($start+$end)/2);
			if(!$score){ $score = $end-$start; }
			push(@{$coord{$chr}},[$pos,$score,$ifile]);
			$sum += $score; 
			$nn++;
			if($chr eq "1"){ die "Err: chr $chr in $filename\n"; }
			if($chr eq "chr20"){ $found=1; }
		}
		$numPeaks[$ifile] = $nn;
		if(!$found){ print "Err: chr20 not found in $filename\n"; }
		if($nn){ $scoreAvg[$ifile] = $sum/$nn; }
	}
	foreach my $chr (sort keys %coord){
		my $ref1 = $coord{$chr};
		@$ref1 = sort {$a->[0]<=>$b->[0]} @$ref1;
		for(my $i=0; $i<@$ref1; $i++){
			my $istart = $i;
			my @positions = ($ref1->[$i]->[0]);
			my @counts;
			for(my $i1=0; $i1<$nfiles; $i1++){ $counts[$i1]=0; }
			$counts[$ref1->[$i]->[2]]=1;
			my $maxScore=$ref1->[$i]->[1]/$scoreAvg[$ref1->[$i]->[2]];
			while($i<@$ref1-1 && $ref1->[$i+1]->[0] - $ref1->[$i]->[0] <= $delta){
				$i++;
				push(@positions,$ref1->[$i]->[0]);
				if(!$counts[$ref1->[$i]->[2]]){
					$counts[$ref1->[$i]->[2]]=1;
				}
				if($maxScore < $ref1->[$i]->[1]/$scoreAvg[$ref1->[$i]->[2]]){
					$maxScore=$ref1->[$i]->[1]/$scoreAvg[$ref1->[$i]->[2]];
				}
			}
			my $median = int(median(\@positions));
			my ($bin,$code)=(1,0);
			my %hashSupport;
			for(my $ifile=0; $ifile<$nfiles; $ifile++){
				$code += $counts[$ifile]*$bin;
				$bin*=2;
				if($counts[$ifile]){
					$hashSupport{$ref->[$ifile]->[1]}++;
				}
			}
			my $Nsupport = keys %hashSupport;
			my $x1 = $Nsupport;
			if($x1>6){ $x1=6; }
			my $score = int(1000*(($x1-1)*100+$maxScore))/1000;
			push(@peaks,[$chr,$median,$score,$Nsupport,$code]);
		}
	}
	@peaks = sort {$b->[2]<=>$a->[2]} @peaks;
	my $nfinal;
	my @freq;
	my @freqTSS;
	my $count=0;
	for(my $i=0; $i<@peaks; $i++){
		my ($chr,$pos,$score,$Nsupport,$code) = @{$peaks[$i]};
		if($i==15000){
			$nfinal=$Nsupport;
			if($nfinal>3){ $nfinal=3; }
		}
		elsif($i>15000 && ($Nsupport<$nfinal || $nfinal==0)){ last; }
		my $refTSS = $hashTSS{$chr};
		my @targets=();
		my @distances=();
		if($refTSS){
			my $NNN = @$refTSS;
			my $iStart = find_in_list($pos-$maxDist,$refTSS, $NNN);
			if($iStart<0){ ++$iStart; }
			#print "$pos\t$refTSS->[$iStart]->[0]\n";
			my $iEnd = $iStart;
			while($iEnd<$NNN && $refTSS->[$iEnd]->[0]<$pos+$maxDist){ ++$iEnd; }
			if($iEnd==$NNN || $refTSS->[$iEnd]->[0]>=$pos+$maxDist){ 
				$iEnd--;
			}
			#print "$iStart $iEnd $pos $refTSS->[$iStart]->[0] $refTSS->[$iEnd]->[0]\n";
			my @score=();
			for(my $it=$iStart; $it<=$iEnd; ++$it){
				my ($TSS,$strand,$symbol) = @{$refTSS->[$it]};
				my $dist = abs($TSS-$pos)/1000;
				if($dist<1){ $dist=1; }
				my $symbolQual = 3;
				if($symbol=~/\d\d\d\d|^C[\dXY]+orf|^FAM\d|^HIST|^MIR\d|^MRP[LS]|^OR\d|^RP[LS]|^SN[OA]R/){ $symbolQual = 1; }
				$score[$it-$iStart] = [$symbolQual/$dist,$it];
			}
			if(@score){
				my %hashDone=();
				my @sorted = sort {$b->[0]<=>$a->[0]} @score;
				my $sign=-1;
				my $itBest = $sorted[0]->[1];
				$hashDone{$refTSS->[$itBest]->[2]} = 1;
				my $dist = $pos - $refTSS->[$itBest]->[0];
				if($dist >= 0){ $sign=1; }
				my $absdist = abs($dist);
				if($absdist < $maxDist){
					push(@targets,$refTSS->[$itBest]->[2]);
					push(@distances,$dist);
					$freqTSS[int($absdist/500)]++;
				}
				for(my $i=1; $i<@sorted; $i++){
					if(keys %hashDone > 2){ last; }
					if($sorted[$i]->[0] >= 0.2*$sorted[0]->[0]){
						my $symbol1 = $refTSS->[$sorted[$i]->[1]]->[2];
						my $dist1 = $pos - $refTSS->[$sorted[$i]->[1]]->[0];
						if($hashDone{$symbol1}){ next; }
						$hashDone{$symbol1} = 1;
						if(abs($dist1) < $maxDist){
							push(@targets,$symbol1);
							push(@distances,$dist1);
						}
					}
				}
			}
		}
		if($chr !~ /^chr/){ next; }
		print OUT "$TF\t$chr\t$pos\t$score\t$Nsupport\t$code";
		if(@targets){
			print OUT "\t".join(',',@targets)."\t".join(',',@distances);
		}
		print OUT "\n";
		if($Nsupport>1){
			my $ifile=0;
			while($code){
				if($code%2){ $freq[$ifile]++; }
				$ifile++;
				$code = int($code/2);
			}
		}
		$count++;
	}
	for(my $ifile=0; $ifile<$nfiles; $ifile++){
		push(@report,[$TF,$ref->[$ifile]->[0],$numPeaks[$ifile],$freq[$ifile],$count,$nfiles,$nfinal]);
		#print "$TF\t$ref->[$ifile]->[0]\t$freq[$ifile]\t$count\t$nfiles\n";
	}
	if(@freqTSS > 15){ splice(@freqTSS,15); }
	#print OUT join("\t",$TF,@freqTSS)."\n";
	push(@report,[$TF,@freqTSS]);
}
for(my $i=0; $i<@report; $i++){
	my @dataReport = @{$report[$i]};
	print OUT join("\t",@dataReport)."\n";
}
close OUT;
exit(0);


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
if(!$n){ return $MISSING; }
my $i = $n/2;
if($i > int($i)){
	$i = int($i);
	$median = $sorted[$i];
}else{
	$median = ($sorted[$i] + $sorted[$i-1])/2;
}
return $median;
}

#*****************
sub find_in_list
#*****************
{
my $start = shift;
my $ref = shift;
my $N = shift;
my $pos = shift;

if(!$pos){ $pos=0; }
if($N<=1){ return 0; }
my $n = $N;
my $j = -1;
my $i = int($n/2);
while(1){
	my $start1=$ref->[$i]->[$pos];
	if($start1==$start){
		return $i;
	}
	elsif($start > $start1){
		if($n-$i<=1){ return $i; }
		$j = $i;
		$i = int(($n+$i)/2);
	}
	else{
		if($i-$j<=1){ if($i){--$i;} return $i; }
		$n = $i;
		$i = int(($j+$n)/2);
	}
}
return $i;
}


