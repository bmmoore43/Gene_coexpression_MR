#!/group/bioinfo/apps/apps/perl-5.20.1/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Math::Round;

###############################################
# Run: perl prefilter_matrix.pl -i [infile] -o [outfile] [options]
#
# prefilter_matrix.pl 
#
#
# Jen Wisecaver
# 20180904, jwisecav@purdue.edu
################################################

# usage
my %opts; getopt('iomn', \%opts );

my $usage = "\n\tUsage Error:
	Run: perl prefilter_matrix.pl -i [infile] -o [outfile] [options]
	
	OPTIONS =
	Files:
	-i: <FULL_FILENAME> INFILE: path input file: matrix of gene counts (Header should start with a 'tab' see example_matrix.txt input file)
	-o: <FULL_FILENAME> OUTFILE: path to outfile: filtered matrix of gene counts 
	
	General:
	-n: <INTEGER> MIN_READS: minimum number of reads mapping to a gene in order for it to be considered 'expresssed' (default = 10)
	-m: <INTEGER> MIN_CONDITIONS: minimum number of conditions in which a gene is expressed in order to retain gene in filtered matrix (default = 3)

	\n\n";


#####################################################
###### Define working files and user variables ######
#####################################################
my $INFILE;
if ($opts{i}) { 
	$INFILE = "$opts{i}";
}
else { die "\n\tInput <i>$usage"; }
my $INPATH = $INFILE;
if ($INPATH =~ /\//){
	$INPATH =~ /\/([^\/]+)$/;
	$INFILE = $1;
}

my $OUTFILE;
if ($opts{o}) { 
	$OUTFILE = "$opts{o}";
}
else { die "\n\tInput <o>$usage"; }
my $OUTPATH = $OUTFILE;
if ($OUTPATH =~ /\//){
	$OUTPATH =~ /\/([^\/]+)$/;
	$OUTFILE = $1;
}

my $MIN_READS = 10; 
if ($opts{n}) { if ($opts{n} =~ m/^\d+$/) { 
	if (int($opts{n}) eq $opts{n}) { $MIN_READS = "$opts{n}"; }} else { die "\n\tInput <n>\n$usage"; } }

my $MIN_CONDITIONS = 3; 
if ($opts{'m'}) { if ($opts{'m'} =~ m/^\d+$/) { 
	if (int($opts{'m'}) eq $opts{'m'}) { $MIN_CONDITIONS = "$opts{'m'}"; }} else { die "\n\tInput <m>\n$usage"; } }


##########################
###### Parse matrix ######
##########################
open (IFIL, '<', $INPATH) or die "Couldn't open file $INPATH: $!\n";
open (OFIL, '>', $OUTPATH) or die "Couldn't write to file $OUTPATH: $!\n";
my $total = 0;
my $saved = 0;
while ( my $line = <IFIL> ) {
	chomp $line;
	if ($line =~ s/^\t//){
		print OFIL "$line\n";
		next;
	}
	
	$total++;
	my @col = split(/\t/, $line);
	my $gene = shift @col;

	my $num = 0;
	my @new_col;
	foreach my $raw_count (@col){
		if ($raw_count >= $MIN_READS){
			$num++;
		}
		my $round_count = round($raw_count);
		push @new_col, $round_count;
		
	}
	
	if ($num >= $MIN_CONDITIONS){
		my $new_line = join("\t", @new_col);
		#print OFIL "$line\n";
		print OFIL "$gene\t$new_line\n";
		$saved++;
	}
}

print "Filtering completed!\n\tOriginal matrix ($INFILE) contained $total genes\n\tFiltered matrix ($OUTFILE) contains $saved genes\n";

close IFIL;
close OFIL;
