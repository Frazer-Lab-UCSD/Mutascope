#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);


my $bed;
my $fasta;
my $pdir;
my $samtools;
my $length=151;

GetOptions("bed=s" => \$bed,
						"fasta=s" => \$fasta,
						"pdir=s" => \$pdir,
						"samtools=s" => \$samtools,
						"length=i" => \$length);


my @rgs;
my %beds;
my %hash;
my @chroms = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM");
my $sample;

open(BED, $bed);
while(my $line = <BED>){
	chomp $line;
	next if($line =~ m/^track/);
	next if($line =~ m/^#/);
	next if($line =~ m/^>/);	
	my ($chr, $bstart, $bend, $name, $score, $strand, $pstart, $pend, @other) = split(/\t/, $line);
	$beds{$chr}{$bstart}{$bend} = $pstart."\t".$pend;
}
close(BED);

my $header = "";

open(SAM, $ARGV[0]);
while(my $line = <SAM>){
	chomp $line;
	if($line =~ m/^@/){
		$header .="$line\n";
	}
	
	if($line =~ m/^\@RG/){
		my @fields = split(/\t/, $line);
		my $rg = $fields[1];
		$rg =~ s/ID://;
		push(@rgs, $rg);
	}
	
	next if($line =~ m/^@/);
	
	my @info = split(/\t/, $line);
	my $flag = $info[1];
	my $tmp_rg;
	foreach (@info){
		if($_ =~ m/^RG:Z:/){
			$tmp_rg = $_;
			$tmp_rg =~ s/^RG:Z://;
		}
	}
	
	my ($samp, $read, $loc) = split(/_/, $tmp_rg);
	$sample = $samp;
	
	if(exists($hash{$tmp_rg})){
		$hash{$tmp_rg} .= "$line\n";
	}else{
		$hash{$tmp_rg} = "$line\n";
	}
}
close(SAM);


my $tmp_counter = 1;
my @vals = keys %hash;
my $len = $#vals;

while($len >= 0){
	my %cur_loc;
	my $cur_lines = "";
	
	open(TMP_OUT, ">$pdir/intermediate/tmp_sam_file_$sample.sam");
	open(BED_OUT, ">$pdir/intermediate/tmp_bed_file_$sample.bed");

	foreach my $chr (@chroms){
		foreach my $start (sort {$a <=> $b} keys %{$beds{$chr}}){
LOOP:			foreach my $stop (sort {$a <=> $b} keys %{$beds{$chr}{$start}}){
				if(exists($hash{"$sample\_R1_$chr:$start-$stop"})){

					if(!(exists($cur_loc{$chr}))){
						$cur_loc{$chr}{$start}{$stop} = 1;
						$cur_lines .= $hash{"$sample\_R1_$chr:$start-$stop"};
						print BED_OUT "$chr\t$start\t$stop\tR1\n";
						delete($hash{"$sample\_R1_$chr:$start-$stop"});
						next LOOP;
					}else{
						next LOOP if(exists($cur_loc{$chr}{$start}{$stop}));
						my $flag = 0;
						foreach my $st (keys %{$cur_loc{$chr}}){
							foreach my $sp (keys %{$cur_loc{$chr}{$st}}){
								$flag = 1 if(($start >= $st and $start <= $sp) or ($stop >= $st and $stop <= $sp));
							}
						}
						if($flag == 0){
							$cur_loc{$chr}{$start}{$stop} = 1;
							$cur_lines .= $hash{"$sample\_R1_$chr:$start-$stop"};
							print BED_OUT "$chr\t$start\t$stop\tR1\n";
							delete($hash{"$sample\_R1_$chr:$start-$stop"});
							next LOOP;
						}
					}
					
				}elsif(exists($hash{"$sample\_R2_$chr:$start-$stop"})){
					if(!(exists($cur_loc{$chr}))){
						$cur_loc{$chr}{$start}{$stop} = 1;
						$cur_lines .= $hash{"$sample\_R2_$chr:$start-$stop"};
						print BED_OUT "$chr\t$start\t$stop\tR2\n";
						delete($hash{"$sample\_R2_$chr:$start-$stop"});
						next LOOP;
					}else{
						next LOOP if(exists($cur_loc{$chr}{$start}{$stop}));
						my $flag = 0;
						foreach my $st (keys %{$cur_loc{$chr}}){
							foreach my $sp (keys %{$cur_loc{$chr}{$st}}){
								$flag = 1 if(($start >= $st and $start <= $sp) or ($stop >= $st and $stop <= $sp));
							}
						}
						if($flag == 0){
							$cur_loc{$chr}{$start}{$stop} = 1;
							$cur_lines .= $hash{"$sample\_R2_$chr:$start-$stop"};
							print BED_OUT "$chr\t$start\t$stop\tR2\n";
							delete($hash{"$sample\_R2_$chr:$start-$stop"});
							next LOOP;
						}
					}
				}
			}
		}
	}
	@vals = keys %hash;
	$len = $#vals;
	
	print TMP_OUT "$header";
	print TMP_OUT "$cur_lines";	
	close(TMP_OUT);
	close(BED_OUT);

	`$samtools view -Sbh $pdir/intermediate/tmp_sam_file_$sample.sam 2> $pdir/intermediate/tmp_stderr_file.txt | $samtools mpileup -Bd 5000000 -Q 0 -O -s -f $fasta -l $pdir/intermediate/tmp_bed_file_$sample.bed - > $pdir/intermediate/tmp_pileup_file_$sample.txt 2> $pdir/intermediate/tmp_stderr_file.txt`;

	`perl $abs_path/convert_pileup_to_xpileup.pl -sample $sample -bed $bed -regions $pdir/intermediate/tmp_bed_file_$sample.bed -length $length -tmp_counter $tmp_counter $pdir/intermediate/tmp_pileup_file_$sample.txt`;

	$tmp_counter++;
	
}

`perl $abs_path/combine_xpileup_files.pl $pdir/intermediate/tmp_pileup_file_$sample\_* > $pdir/intermediate/$sample.xpileup`;

`rm $pdir/intermediate/tmp_pileup_file_$sample.txt`;
`rm $pdir/intermediate/tmp_bed_file_$sample.bed`;
`rm $pdir/intermediate/tmp_pileup_file_$sample\_*`;
`rm $pdir/intermediate/tmp_sam_file_$sample.sam`;
`rm $pdir/intermediate/tmp_stderr_file.txt`;




