#! /usr/bin/perl -w

use warnings;
use strict;

my %hash;

open(FILE, $ARGV[0]);
while(my $line = <FILE>){
	chomp $line;
	next if($line =~ m/^track/);
	next if($line =~ m/^#/);
	next if($line =~ m/^>/);	
	my ($chr, $start, $stop, @other) = split(/\t/, $line);
	print "$chr:$start-$stop\n";
}

