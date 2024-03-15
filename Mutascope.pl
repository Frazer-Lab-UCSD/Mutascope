#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);

my $Mutascope_usage = "\nMutascope is a software package that is designed to analyze and process PCR amplicon deep sequencing data.

Mutascope COMMANDS:
\trunPipeline\tRun the full Mutascope pipeline (Recommended)
\tmakeBlackList\tPrepared the list of excluded mutations based on homologous sequence
\trunBWA\t\tAlignes reads with BWA's Smith-Waterman aligner
\trefinement\tRefigns the aligned reads by removing low quality mapped reads, multi-mapped reads, and more
\tgroupRealign\tRealigns a Tumor-Normal pair together using GATK to improve SNP and Indel calling
\txpileup\t\tGenerates the an extended pileup (xpileup) file used to call variants and calculate error rates
\tcalcErrorRates\tCalculates the error rates given the position in the read, read orientation, and mutation type
\tcallSomatic\tCalls somatic variants and applies a set of standard filters\n";


die "$Mutascope_usage\n" if(!$ARGV[0]);

my $command = shift(@ARGV);
my $exit_status;

if($command eq "runPipeline"){
	$exit_status = system("perl $abs_path/scripts/runPipeline.pl @ARGV");
	die "\n" if($exit_status != 0);
}elsif($command eq "makeBlackList"){
	$exit_status = system("perl $abs_path/scripts/makeBlackList.pl @ARGV");
	die "\n" if($exit_status != 0);
}elsif($command eq "runBWA"){
	$exit_status = system("perl $abs_path/scripts/runBWA.pl @ARGV");
	die "\n" if($exit_status != 0);
}elsif($command eq "refinement"){
	$exit_status = system("perl $abs_path/scripts/refinement.pl @ARGV");
	die "\n" if($exit_status != 0);
}elsif($command eq "groupRealign"){
	$exit_status = system("perl $abs_path/scripts/groupRealign.pl @ARGV");
	die "\n" if($exit_status != 0);
}elsif($command eq "xpileup"){
	$exit_status = system("perl $abs_path/scripts/xpileup.pl @ARGV");
	die "\n" if($exit_status != 0);
}elsif($command eq "calcErrorRates"){
	$exit_status = system("perl $abs_path/scripts/calcErrorRates.pl @ARGV");
	die "\n" if($exit_status != 0);
}elsif($command eq "callSomatic"){
	$exit_status = system("perl $abs_path/scripts/callSomatic.pl @ARGV");
	die "\n" if($exit_status != 0);
}else{
	die "\n$Mutascope_usage\nUnknown command \"$command\"\n\n";
}






