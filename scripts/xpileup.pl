#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);

sub check_for_samtools();
sub check_outDir($);
sub check_bed($);

my $usage = "\nUse:\n\tMutascope xpileup [-h] [-bed FILE] [-fasta FILE] [-pdir PATH] [-sample name] [-samtools_path FILE]

Description:
xpileup generates the extended pileup file (xpileup) required by Mutascope callSomatic (in gzip format). The xpileup contains information about read groups and nucleotide positions in the reads. The program requires SAMTools to be installed.

Example:
\tMutascope xpileup -bed target_primers.bed -fasta hg19.fa -pdir /home/projects/ -sample SAMPLE

Options:
\t-bed\t\tBED\tBED format file containing target locations with primer locations (See BED format in the Manual.pdf for more detail) [REQUIRED]
\t-fasta\t\tFA\tFASTA file of reference genome [REQUIRED]
\t-pdir\t\tPATH\tPath to the PROJECT directory [REQUIRED]
\t-sample\t\tSTRING\tName of the sample (only use letters or numbers! No symbols!) [REQUIRED]
\t-samtools_path\tFILE\tSAMTools executable (full path) [Optional]
\t-length\t\tINT\tLength of the sequencing reads. If LENGTH varies than input the length of the longest read (Default 151) [Optional]
\t-h\t\tFLAG to output help message [Optional]
\n";


my $bed;
my $fasta;
my $samtools_path;
my $pdir;
my $sample;
my $h;
my $length=151;

GetOptions("bed=s" => \$bed,
						"fasta=s" => \$fasta,
						"samtools_path=s" => \$samtools_path,
						"pdir=s" => \$pdir,
						"sample=s" => \$sample,
						"length=i" => \$length,
						"h" => \$h);


die "$usage" if($h);
die "\n$usage\nERROR: Must provide BED file!\n\n" if(!$bed);
die "\n$usage\nERROR: Bed file \"$bed\" doesn't exist!\n\n" if(!(-e $bed));
die "\n$usage\nERROR: Must provide FASTA file!\n\n" if(!$fasta);
die "\n$usage\nERROR: Fasta file \"$fasta\" doesn't exist!\n\n" if(!(-e $fasta));
die "\n$usage\nERROR: Must provide the full path to the PROJECT Directory!\n\n" if(!$pdir);
die "\n$usage\nERROR: The Directory \"$pdir\" doesn't exists!\n\n" if(!(-d $pdir));
die "\n$usage\nERROR: Please Specify the name of the sample. This should be consistent throughout all Mutascope modules.\n\n" if(!$sample);
die "\n$usage\nERROR: The sample name can NOT have a '_' in the name.  Please use a different sample name other than \"$sample\".\n\n"if($sample =~ m/_/);

my $bam;
if(-e "$pdir/intermediate/$sample\_merged_unClipped_primersClipped_realigned.bam"){
	$bam = "$pdir/intermediate/$sample\_merged_unClipped_primersClipped_realigned.bam";
}elsif(-e "$pdir/intermediate/$sample\_merged_unClipped_primersClipped.bam"){
	$bam = "$pdir/intermediate/$sample\_merged_unClipped_primersClipped.bam";
}else{
	die "\n$usage\nERROR: The file $pdir/intermediate/$sample\_merged_unClipped_primersClipped_realigned.bam   OR   $pdir/intermediate/$sample\_merged_unClipped_primersClipped.bam  do not exists!\nCheck to make sure you put in the correct sample name ($sample) and PROJECT directory ($pdir).\n\nThis module assumes you previously ran groupRealign or refinement modules.\n\n";
}

check_outDir($pdir);
check_bed($bed);

if(!$samtools_path){
	$samtools_path = check_for_samtools();
}else{
	die "\n$usage\nERROR: SAMTools: \"$samtools_path\" doesn't exists!\n\n" if(!(-e $samtools_path));
}
print STDERR "Using SAMTools: $samtools_path\n";





`$samtools_path view -h $bam | perl $abs_path/make_xpileup_file.pl -length $length -bed $bed -fasta $fasta -pdir $pdir -samtools $samtools_path - `;
`gzip $pdir\/intermediate/$sample.xpileup`;





sub check_bed($){
	my ($bed) = @_;
	open(BED, $bed);
	while(my $line = <BED>){
		chomp $line;
		next if($line =~ m/^track/);
		next if($line =~ m/^#/);
		next if($line =~ m/^>/);
		
		my ($chr, $start, $stop, $name, $score, $strand, $pstart, $pstop, @other) = split(/\t/, $line);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if(!$pstart);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if(!$pstop);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if(!$start);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if(!$stop);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if($start > $pstart);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if($stop < $pstop);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if($pstart > $pstop);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if($start > $stop);
		
	}
	close(BED);
}

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








