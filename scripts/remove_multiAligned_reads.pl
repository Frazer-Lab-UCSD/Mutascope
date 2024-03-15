#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;

my @reads;
my $prev_rname = "tmp";
my $total = 0;
my $mapped = 0;
my $count = 0;

my $read;

GetOptions("read=s" => \$read);


open(SAM, $ARGV[0]);
while(my $line = <SAM>){
	chomp $line;
	print "$line\n" if($line =~ m/^@/);
	next if($line =~ m/^@/);
	
	$total++;
	
	my ($rname, $flag, @other) = split(/\t/, $line);
	
	if($flag == 4){
		next;
	}
	
	if($prev_rname eq "tmp"){
		$prev_rname = $rname;
		push(@reads, $line);
	}else{
		if($rname eq $prev_rname){
			push(@reads, $line);
			$total--;
		}else{
			$mapped++;
			if($#reads == 0){
				$prev_rname = $rname;
				print "$reads[0]\n";
				$count++;
				@reads = ();
				push(@reads, $line);
			}else{
				$prev_rname = $rname;
				@reads = ();
				push(@reads, $line);
			}
		}
	}
}

$mapped++;
if($#reads == 0){
	print "$reads[0]\n";
	$count++;
}

open(OUT, ">>$ARGV[1]");

print OUT "$read\_sequenced_reads\t$total\n";
print OUT "$read\_mapped_reads\t$mapped\n";
print OUT "$read\_uniquely_aligned_reads\t$count\n";

close (OUT);
