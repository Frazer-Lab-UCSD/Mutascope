#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;

my %bed_start;
my %bed_stop;
my %strands;

my $sample;
my $read;
my $dist = 2;
my $strand_specific = "FALSE";

GetOptions("sample=s" => \$sample,
						"read=s" => \$read,
						"dist=i" => \$dist,
						"strand_specific=s" => \$strand_specific);

my @rgs;

open(BED, $ARGV[0]);
while(my $line = <BED>){
	chomp $line;
	next if($line =~ m/^track/);
	next if($line =~ m/^#/);
	next if($line =~ m/^>/);
	my ($chr, $start, $stop, $gene, $score, $strand, @other) = split(/\t/, $line);
	$bed_start{$chr}{$start} = $stop;
	$bed_stop{$chr}{$stop} = $start;
	$strands{$chr}{$start}{$stop} = $strand;
	push(@rgs, "$sample\_$read\_$chr:$start-$stop");
}

my $rg_flag = 0;

my $count = 0;
my $total = 0;
my $indel_reads = 0;
my $soft_clipped_reads = 0;


open(SAM, $ARGV[1]);
READLINE: while(my $line = <SAM>){
	chomp $line;

	if($line =~ m/^@/){
		print "$line\n";
		next READLINE;
	}
	if($rg_flag == 0){
		$rg_flag = 1;
		foreach my $rgs_name (@rgs){
			print "\@RG\tID:$rgs_name\tPL:illumina\tPU:\tLB:$sample.LX1\tDS:$rgs_name\tSM:$sample\n";
		}
	}

	my @fields = split(/\t/, $line);
	my $flag = $fields[1];
	my $chr = $fields[2];
	my $start = $fields[3];
	my $cigar = $fields[5];
	my $rg = $sample."_".$read;
	
	$total++;

	my $cur_start = $start;
	my $cur_end = $start-1;

	if($cigar =~ m/^\d+S/){
		my ($len, @tmps) = split(/S/, $cigar);
		$cur_start -= $len;
	}
	if($cigar =~ m/(H|S|I|D)/){
		my @cigar_counts = split(/\D+/, $cigar);
		my @cigar_types = split(/\d{1,3}/, $cigar);
		shift @cigar_types;
		
		for(my $i = 0; $i <= $#cigar_types; $i++){
			my $type = $cigar_types[$i];
			my $count = $cigar_counts[$i];
			if($type eq "M" or $type eq "D"){
				$cur_end += $count;
			}
		}
		if($cigar_types[$#cigar_types] eq "S"){
			$cur_end += $cigar_counts[$#cigar_types];
		}
	}else{
		$cigar =~ s/M//;
		$cur_end += $cigar;
	}
	
	if(exists($bed_start{$chr}{$cur_start}) and exists($bed_stop{$chr}{$cur_end})){
		if($strand_specific eq "TRUE"){
			if(($read eq "R1" and $strands{$chr}{$cur_start}{$bed_start{$chr}{$cur_start}} eq "+" and $flag == 0) or ($read eq "R1" and $strands{$chr}{$bed_stop{$chr}{$cur_end}}{$cur_end} eq "-" and $flag == 16)){
				if($read eq "R1" and $strands{$chr}{$cur_start}{$bed_start{$chr}{$cur_start}} eq "+" and $flag == 0){
					$rg .= "_$chr:$cur_start-$bed_start{$chr}{$cur_start}";				
				}else{
					$rg .= "_$chr:$bed_stop{$chr}{$cur_end}-$cur_end";				
				}
			}elsif(($read eq "R2" and $strands{$chr}{$cur_start}{$bed_start{$chr}{$cur_start}} eq "+" and $flag == 16) or ($read eq "R2" and $strands{$chr}{$bed_stop{$chr}{$cur_end}}{$cur_end} eq "-" and $flag == 0)){
				if($read eq "R2" and $strands{$chr}{$cur_start}{$bed_start{$chr}{$cur_start}} eq "+" and $flag == 16){
					$rg .= "_$chr:$cur_start-$bed_start{$chr}{$cur_start}";				
				}else{
					$rg .= "_$chr:$bed_stop{$chr}{$cur_end}-$cur_end";				
				}
			}else{
				next READLINE;
			}
		}else{
			my $rand = rand();
			if($rand <= 0.5){
				$rg .= "_$chr:$cur_start-$bed_start{$chr}{$cur_start}";				
			}else{
				$rg .= "_$chr:$bed_stop{$chr}{$cur_end}-$cur_end";				
			}
		}
	}elsif(exists($bed_start{$chr}{$cur_start})){
		if($strand_specific eq "TRUE"){
			if($read eq "R1"){
				next READLINE if($strands{$chr}{$cur_start}{$bed_start{$chr}{$cur_start}} eq "+" and $flag == 16); 
				next READLINE if($strands{$chr}{$cur_start}{$bed_start{$chr}{$cur_start}} eq "-" and $flag == 0); 
			}else{
				next READLINE if($strands{$chr}{$cur_start}{$bed_start{$chr}{$cur_start}} eq "+" and $flag == 0); 
				next READLINE if($strands{$chr}{$cur_start}{$bed_start{$chr}{$cur_start}} eq "-" and $flag == 16); 
			}
		}
		$rg .= "_$chr:$cur_start-$bed_start{$chr}{$cur_start}";			
	}elsif(exists($bed_stop{$chr}{$cur_end})){
		if($strand_specific eq "TRUE"){
			if($read eq "R1"){
				next READLINE if($strands{$chr}{$bed_stop{$chr}{$cur_end}}{$cur_end} eq "+" and $flag == 16); 
				next READLINE if($strands{$chr}{$bed_stop{$chr}{$cur_end}}{$cur_end} eq "-" and $flag == 0); 
			}else{
				next READLINE if($strands{$chr}{$bed_stop{$chr}{$cur_end}}{$cur_end} eq "+" and $flag == 0); 
				next READLINE if($strands{$chr}{$bed_stop{$chr}{$cur_end}}{$cur_end} eq "-" and $flag == 16); 
			}
		}
		$rg .= "_$chr:$bed_stop{$chr}{$cur_end}-$cur_end";				
	}else{
	
		my $min_dist_start = $dist+1;
		my $min_bed_start = "tmp";
		my $start_inc_flag = 0;
		foreach my $st (keys %{$bed_start{$chr}}){
			if($min_dist_start > abs($st - $cur_start)){
				if($start_inc_flag == 0){
					$min_bed_start = $st;
					$min_dist_start = abs($st - $cur_start);
					$start_inc_flag = 1 if($st <= $cur_start and $bed_start{$chr}{$st} >= $cur_end);
				}elsif($st <= $cur_start and $bed_start{$chr}{$st} >= $cur_end){
					$min_bed_start = $st;
						$min_dist_start = abs($st - $cur_start);							
				}
			}
		}
		
		my $min_dist_end = $dist+1;
		my $min_bed_end = "tmp";
		my $end_inc_flag = 0;
		foreach my $sp (keys %{$bed_stop{$chr}}){
			if($min_dist_end > abs($sp - $cur_end)){
				if($end_inc_flag == 0){
					$min_bed_end = $sp;
					$min_dist_end = abs($sp - $cur_end);
					$end_inc_flag = 1 if($bed_stop{$chr}{$sp} <= $cur_start and $sp >= $cur_end);
				}elsif($bed_stop{$chr}{$sp} <= $cur_start and $sp >= $cur_end){
					$min_bed_end = $sp;
					$min_dist_end = abs($sp - $cur_end);							
				}
			}
		}
		
		next READLINE if($min_dist_start > $dist and $min_dist_end > $dist);
		
		if($start_inc_flag == 1 and $end_inc_flag == 1){
			if($min_dist_start <= $min_dist_end){
				if($strand_specific eq "TRUE"){
					if($read eq "R1"){
						next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "+" and $flag == 16); 
						next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "-" and $flag == 0); 
					}else{
						next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "+" and $flag == 0); 
						next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "-" and $flag == 16); 
					}
				}
				$rg .= "_$chr:$min_bed_start-$bed_start{$chr}{$min_bed_start}";
			}else{
				if($strand_specific eq "TRUE"){
					if($read eq "R1"){
						next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "+" and $flag == 16); 
						next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "-" and $flag == 0); 
					}else{
						next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "+" and $flag == 0); 
						next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "-" and $flag == 16); 
					}
				}
				$rg .= "_$chr:$bed_stop{$chr}{$min_bed_end}-$min_bed_end";						
			}
		}elsif($start_inc_flag == 1){
			if($strand_specific eq "TRUE"){
				if($read eq "R1"){
					next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "+" and $flag == 16); 
					next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "-" and $flag == 0); 
				}else{
					next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "+" and $flag == 0); 
					next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "-" and $flag == 16); 
				}
			}
			$rg .= "_$chr:$min_bed_start-$bed_start{$chr}{$min_bed_start}";					
		}elsif($end_inc_flag == 1){
			if($strand_specific eq "TRUE"){
				if($read eq "R1"){
					next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "+" and $flag == 16); 
					next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "-" and $flag == 0); 
				}else{
					next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "+" and $flag == 0); 
					next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "-" and $flag == 16); 
				}
			}
			$rg .= "_$chr:$bed_stop{$chr}{$min_bed_end}-$min_bed_end";											
		}else{
			if($min_dist_start <= $min_dist_end){
				if($strand_specific eq "TRUE"){
					if($read eq "R1"){
						next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "+" and $flag == 16); 
						next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "-" and $flag == 0); 
					}else{
						next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "+" and $flag == 0); 
						next READLINE if($strands{$chr}{$min_bed_start}{$bed_start{$chr}{$min_bed_start}} eq "-" and $flag == 16); 
					}
				}			
				$rg .= "_$chr:$min_bed_start-$bed_start{$chr}{$min_bed_start}";
			}else{
				if($strand_specific eq "TRUE"){
					if($read eq "R1"){
						next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "+" and $flag == 16); 
						next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "-" and $flag == 0); 
					}else{
						next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "+" and $flag == 0); 
						next READLINE if($strands{$chr}{$bed_stop{$chr}{$min_bed_end}}{$min_bed_end} eq "-" and $flag == 16); 
					}
				}
				$rg .= "_$chr:$bed_stop{$chr}{$min_bed_end}-$min_bed_end";						
			}
		}
	}
	
	print "$line\tRG:Z:$rg\n";
	$count++;
	$indel_reads++ if($fields[5] =~ m/I|D/);
	$soft_clipped_reads++ if($fields[5] =~ m/S/);
}

open(OUT, ">>$ARGV[2]");
print OUT "$read\_reads_with_expected_alignment\t$count\n";
print OUT "\t$read\_expected_alignment_reads_with_an_indel\t$indel_reads\n";
print OUT "\t$read\_expected_alignment_reads_with_soft_clipping\t$soft_clipped_reads\n";
close(OUT);
