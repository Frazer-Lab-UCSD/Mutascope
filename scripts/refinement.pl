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

my $usage = "\nUse:\n\tMutascope refinement [-h] [-bed BED] [-pdir PATH] [-sample NAME] [-samtools_path FILE] [-dist INT] [-strand_specific BOOL] 

Description:
\trefinement removes low quality reads and edits remaining alignments.  

Example:
\tMutascope refinement -bed BED_NAME.bed -pdir /home/projects/sample/ -sample SAMPLE -strand_specific TRUE -r1_strand positive -r2_strand negative

Options:
\t-bed\t\tBED\tBED format file containing target locations with primer locations (See BED format in the Manual.pdf for more detail) [REQUIRED]
\t-pdir\t\tPATH\tPath to the PROJECT directory [REQUIRED]
\t-sample\t\tSTRING\tName of the sample (only use letters or numbers! No symbols!). This should be consistent throughout all Mutascope modules [REQUIRED]
\t-samtools_path\tFILE\tSAMTools executable (full path) [Optional]
\t-dist\t\tINT\tRead vs Amplicon start or stop alignment offset, in nucleotides (default 2) [Optional]
\t-strand_specific\tTRUE/FALSE\tReads are strand specific e.g. all read1 align on the same strand. (default FALSE) [Optional]
\t-h\t\tFLAG to output help message [Optional]
\n";



my $bed;
my $samtools_path;
my $pdir;
my $sample;
my $n=3;
my $c=5;
my $g=1;
my $e=6;
my $h;
my $dist = 2;
my $strand_specific = "FALSE";


GetOptions("bed=s" => \$bed,
						"samtools_path=s" => \$samtools_path,
						"pdir=s" => \$pdir,
						"sample=s" => \$sample,
						"n=i" => \$n,
						"c=i" => \$c,
						"g=i" => \$g,
						"e=i" => \$e,
						"dist=i" => \$dist,
						"strand_specific=s" => \$strand_specific,
						"h" => \$h);

die "$usage" if($h);
die "\n$usage\nMust provide BED file!\n\n" if(!$bed);
die "\n$usage\nBed file \"$bed\" doesn't exist!\n\n" if(!(-e $bed));
die "\n$usage\nMust provide the full path to the Output Directory!\n\n" if(!$pdir);
die "\n$usage\nThe Directory \"$pdir\" doesn't exists!\n\n" if(!(-d $pdir));
die "\n$usage\nPlease Specify the name of the sample. This should be consistent throughout all Mutascope modules.\n\n" if(!$sample);
die "\n$usage\nThe sample name can NOT have a '_' in the name.  Please use a different sample name other than \"$sample\".\n\n"if($sample =~ m/_/);

check_outDir($pdir);
check_bed($bed);

if(!$samtools_path){
	$samtools_path = check_for_samtools();
}else{
	die "\n$usage\nSAMTools: \"$samtools_path\" doesn't exists!\n\n" if(!(-e $samtools_path));
}
print STDERR "Using SAMTools: $samtools_path\n";

die "\n$usage\nThe file $pdir\/intermediate/$sample\_read1_bwaSWAln.bam doesn't exists!\nCheck to make sure you put in the correct sample name ($sample).\n\nThis module assumes you ran Mutascope runBWA. If you did not you need to first run 'runBWA' before running this module\n\n" if(!(-e "$pdir\/intermediate/$sample\_read1_bwaSWAln.bam"));
die "\n$usage\nThe file $pdir\/intermediate/$sample\_read1_bwaSWAln.bam doesn't exists!\nCheck to make sure you put in the correct sample name ($sample).\n\nThis module assumes you ran Mutascope runBWA. If you did not you need to first run 'runBWA' before running this module\n\n" if(!(-e "$pdir\/intermediate/$sample\_read2_bwaSWAln.bam"));

open(OUT, ">$pdir\/quality/$sample\_alignment_stats.txt");
print OUT "";
close(OUT);

`$samtools_path view -h $pdir\/intermediate/$sample\_read1_bwaSWAln.bam | perl $abs_path/remove_multiAligned_reads.pl -read R1 - $pdir\/quality/$sample\_alignment_stats.txt | perl $abs_path/filter_SAM_on_SW_score.pl -read R1 -n $n -c $c -g $g -e $e - $pdir\/quality/$sample\_alignment_stats.txt | perl $abs_path/add_readGroup_to_SAM_and_filter.pl -sample $sample -read R1 -dist $dist -strand_specific $strand_specific $bed - $pdir\/quality/$sample\_alignment_stats.txt | $samtools_path view -Sbh - | $samtools_path sort - $pdir\/intermediate/$sample\_read1_mapped_filteredSW_rgAddedFilt_sorted`;

`$samtools_path view -h $pdir\/intermediate/$sample\_read2_bwaSWAln.bam | perl $abs_path/remove_multiAligned_reads.pl -read R2 - $pdir\/quality/$sample\_alignment_stats.txt | perl $abs_path/filter_SAM_on_SW_score.pl -read R2 -n $n -c $c -g $g -e $e - $pdir\/quality/$sample\_alignment_stats.txt | perl $abs_path/add_readGroup_to_SAM_and_filter.pl -sample $sample -read R2 -dist $dist -strand_specific $strand_specific $bed - $pdir\/quality/$sample\_alignment_stats.txt | $samtools_path view -Sbh - | $samtools_path sort - $pdir\/intermediate/$sample\_read2_mapped_filteredSW_rgAddedFilt_sorted`;

`mkdir $pdir/intermediate/$$/`;
`java -Xmx2048M -Djava.io.tmpdir=$pdir/intermediate/$$/ -jar $abs_path/MergeSamFiles.jar I=$pdir\/intermediate/$sample\_read1_mapped_filteredSW_rgAddedFilt_sorted.bam I=$pdir\/intermediate/$sample\_read2_mapped_filteredSW_rgAddedFilt_sorted.bam O=$pdir\/intermediate/$sample\_merged.bam AS=TRUE USE_THREADING=TRUE`;
`rm -rf $pdir/intermediate/$$`;

`$samtools_path view -h $pdir\/intermediate/$sample\_merged.bam | perl $abs_path/remove_softClipping.pl - | $samtools_path view -Shb - | $samtools_path sort - $pdir\/intermediate/$sample\_merged_unClipped`;

`$samtools_path view -h $pdir\/intermediate/$sample\_merged_unClipped.bam | perl $abs_path/clip_primers_from_bam.pl $bed - | $samtools_path view -Sbh - | $samtools_path sort - $pdir\/intermediate/$sample\_merged_unClipped_primersClipped`;

`$samtools_path index $pdir\/intermediate/$sample\_merged_unClipped_primersClipped.bam`;

`$samtools_path view -h $pdir\/intermediate/$sample\_merged_unClipped_primersClipped.bam | perl $abs_path/calculate_bedReadCounts.pl $bed - > $pdir\/quality/$sample\_readsPerAmplicon.txt`;

`R --slave --args $pdir\/quality/$sample\_readsPerAmplicon.txt $pdir\/quality/$sample\_readsPerAmplicon_meanSD.txt $pdir\/quality/$sample\_readsPerAmplicon_cov2X.txt $pdir\/quality/$sample\_readsPerAmplicon_sensitivity.txt $pdir\/quality/$sample\_readsPerAmplicon_cumulativeCoverageDist.pdf < $abs_path/read_cov_rcode.r`;


`rm $pdir\/intermediate/$sample\_read1_mapped_filteredSW_rgAddedFilt_sorted.bam`; 
`rm $pdir\/intermediate/$sample\_read2_mapped_filteredSW_rgAddedFilt_sorted.bam`; 
`rm $pdir\/intermediate/$sample\_merged_unClipped.bam`;
`rm $pdir\/intermediate/$sample\_merged.bam`;


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








