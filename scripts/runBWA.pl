#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);

sub check_for_bwa();
sub check_for_samtools();
sub check_bwa_index($);
sub check_outDir($);

my $usage = "\nUse:\n\tMutascope runBWA [-h] [-fasta FA] [-pdir PATH] [-sample NAME] [-samtools_path FILE] [-bwa_path FILE] [-r1_fastq FILE] [-r2_fastq FILE] [-t INT]

Description:
runBWA takes in a set of FASTQ files (r1_fastq and r2_fastq) and aligns the reads to the genome using BWA-SW. It requires BWA and SAMTools to be installed.

Example:
\tMutascope runBWA -fasta hg19.fa -pdir /home/projects/ -sample SAMPLE -r1_fastq SAMPLE_r1.fastq -r2_fastq SAMPLE_r2.fastq -t 10

Options:
\t-fasta\t\tFA\tFASTA file of reference genome [REQUIRED]
\t-pdir\t\tPATH\tPath to the PROJECT directory [REQUIRED]
\t-sample\t\tSTRING\tName of the sample (only use letters or numbers! No symbols!). This should be consistent throughout all Mutascope modules [REQUIRED]
\t-r1_fastq\tFASTQ\tRead1 FASTQ file [REQUIRED]
\t-r2_fastq\tFASTQ\tRead2 FASTQ file [REQUIRED]
\t-samtools_path\tFILE\tSAMTools executable (full path) [Optional]
\t-bwa_path\tFILE\tBWA executable (full path) [Optional]
\t-t\t\tINT\tNumber of threads for BWA (default 1) [Optional]
\t-h\t\tFLAG to output help message [Optional]
\n";



my $fasta;
my $bwa_path;
my $samtools_path;
my $r1_fastq;
my $r2_fastq;
my $t = 1;
my $z = 5;
my $s = 5;
my $pdir;
my $sample;
my $h;

GetOptions("fasta=s" => \$fasta,
						"bwa_path=s" => \$bwa_path,
						"samtools_path=s" => \$samtools_path,
						"r1_fastq=s" => \$r1_fastq,
						"r2_fastq=s" => \$r2_fastq,
						"t=i" => \$t,
						"z=i" => \$z,
						"s=i" => \$s,
						"pdir=s" => \$pdir,
						"h" => \$h,
						"sample=s" => \$sample);


die "$usage" if($h);
die "\n$usage\nERROR: Must provide FASTA file!\n\n" if(!$fasta);
die "\n$usage\nERROR: Fasta file \"$fasta\" doesn't exist!\n\n" if(!(-e $fasta));
die "\n$usage\nERROR: FASTQ file \"$r1_fastq\" for r1_fastq doesn't exists!\n\n" if(!$r1_fastq or !(-e $r1_fastq));
die "\n$usage\nERROR: FASTQ file \"$r2_fastq\" for r2_fastq doesn't exists!\n\n" if(!$r2_fastq or !(-e $r2_fastq));
die "\n$usage\nERROR: Must provide the full path to the Output Directory!\n\n" if(!$pdir);
die "\n$usage\nERROR: The Directory \"$pdir\" doesn't exists!\n\n" if(!(-d $pdir));
die "\n$usage\nERROR: Please Specify the name of the sample. This should be consistent throughout all Mutascope modules.\n\n" if(!$sample);
die "\n$usage\nERROR: The sample name can NOT have a '_' in the name.  Please use a different sample name other than \"$sample\".\n\n"if($sample =~ m/_/);

if(!$bwa_path){
	$bwa_path = check_for_bwa();
}else{
	die "\n$usage\nERROR: BWA: \"$bwa_path\" doesn't exists!\n\n" if(!(-e $bwa_path));
}
print STDERR "Using bwa: $bwa_path\n\n";

if(!$samtools_path){
	$samtools_path = check_for_samtools();
}else{
	die "\n$usage\nERROR: SAMTools: \"$samtools_path\" doesn't exists!\n\n" if(!(-e $samtools_path));
}
print STDERR "Using SAMTools: $samtools_path\n\n";

check_bwa_index($fasta);
check_outDir($pdir);


`$bwa_path bwasw -z $z -t $t -s $s $fasta $r1_fastq | $samtools_path view -Shb - > $pdir\/intermediate/$sample\_read1_bwaSWAln.bam`;
`$bwa_path bwasw -z $z -t $t -s $s $fasta $r2_fastq | $samtools_path view -Shb - > $pdir\/intermediate/$sample\_read2_bwaSWAln.bam`;





sub check_outDir($){
	my ($od) = @_;
	if(!(-d "$od/results")){
		`mkdir $od/results`;
		print STDERR "Making \"$od/results\" directory!\n";
	}
	if(!(-d "$od/intermediate")){
		`mkdir $od/intermediate`;
		print STDERR "Making \"$od/intermediate\" directory!\n";
	}
	if(!(-d "$od/quality")){
		`mkdir $od/quality`;
		print STDERR "Making \"$od/quality\" directory!\n";
	}
}

sub check_for_samtools(){
	my $samtools_results = `which samtools`;
	chomp $samtools_results;
	if($samtools_results eq ""){
		die "\nERROR: Can not find SAMTools. Please install SAMTools (http://samtools.sourceforge.net/).\n   If installed either specify the path (-samtools_path) or create a link to your local path (Example: /usr/bin/).\n"; 
	}else{
		return $samtools_results;
	}
}


sub check_for_bwa(){
	my $bwa_result = `which bwa`;
	chomp $bwa_result;
	if($bwa_result eq ""){
		die "\nERROR: Can not find BWA. Please install BWA (http://bio-bwa.sourceforge.net/).\n   If installed either specify the path (-bwa_path) or create a link to your local path (Example: /usr/bin/).\n"; 
	}else{
		return $bwa_result;
	}
}

sub check_bwa_index($){
	my ($fasta) = @_;
	my $index_flag = 0;
	my $missing = "";
	if(!(-e "$fasta.amb")){
		$index_flag = 1;
		$missing = "$fasta.amb";
	}elsif(!(-e "$fasta.ann")){
		$index_flag = 1;
		$missing = "$fasta.ann";
	}elsif(!(-e "$fasta.bwt")){
		$index_flag = 1;
		$missing = "$fasta.bwt";
	}elsif(!(-e "$fasta.pac")){
		$index_flag = 1;
		$missing = "$fasta.pac";
	}elsif(!(-e "$fasta.sa")){
		$index_flag = 1;
		$missing = "$fasta.sa";
	}
	
	if($index_flag == 1){
		die "\nERROR:\tThe fasta file \"$fasta\" is not indexed properly by BWA. You are missing the \"$missing\" file\n\nPlease index the fasta file with bwa! Run: bwa index -a bwtsw $fasta\n\n";
		#`bwa index -a bwtsw $fasta`;
		
	}
}

