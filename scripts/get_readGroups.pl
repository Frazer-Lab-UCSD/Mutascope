#! /usr/bin/perl -w

use warnings;
use strict;

open(FILE, $ARGV[0]);
while(my $line = <FILE>){
	chomp $line;
	next if(!($line =~ m/^\@RG/));
	my ($rg, $ID, @other) = split(/\t/, $line);
	$ID =~ s/^ID://;
	print "$ID\n";
}
