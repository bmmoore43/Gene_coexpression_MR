#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);

###############################################
# Run: perl calc_mutual_rank.pl -i [infile] -o [outfile] [options]
#
# calc_mutual_rank.pl
#
#
# Jen Wisecaver
# 20180904, jwisecav@purdue.edu
################################################

# usage
my %opts; getopt('iocw', \%opts );

my $usage = "\n\tUsage Error:
	Run: perl calc_mutual_rank.pl -i [infile] -o [outfile] [options]
	
	OPTIONS =
	Files:
	-i: <DIRECTORY> INDIR: path to directory containing PCCs
	-o: <FILENAME> OUTFILE: path to output file
	
	General:
	-c: <NUMBER> PCC_THRESHOLD: retained edges must have a PCC >= this value (default = 0.3)
	-w: <NUMBER> WEIGHT_THRESHOLD: retained edges must have a weight (transformed MR) >= this value (default = 0.01)

	\n\n";


#####################################################
###### Define working files and user variables ######
#####################################################
my $INPATH;
if ($opts{i}) { 
	$INPATH = "$opts{i}";
}
else { die "\n\tInput <i>$usage"; }

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

my $PCC_THRESHOLD = 0.3; 
if ($opts{c}) { if ( looks_like_number($opts{c}) ) { $PCC_THRESHOLD = "$opts{c}"; } else { die "\n\tInput <c>\n$usage"; } }

my $WEIGHT_THRESHOLD = 0.01; 
if ($opts{w}) { if ( looks_like_number($opts{w}) ) { $WEIGHT_THRESHOLD = "$opts{w}"; } else { die "\n\tInput <w>\n$usage"; } }

###########################
###### Parse PCC Dir ######
###########################

my @files = glob("$INPATH/*");
my %bighash;

foreach my $infile (@files) {
	$infile =~ /\/*([^\/]+)$/;
	my $gene1 = $1;
	print "READING: $gene1\n";

	open (IFIL, '<', $infile) or die "Couldn't open file $infile: $!\n";
	while ( my $line = <IFIL> ) {
		chomp $line;
		my ($gene2, $cor) = split(/\t/, $line);
		next if ($gene2 eq $gene1);
		next if (exists $bighash{"$gene2\t$gene1"});
	
		if ($cor > $PCC_THRESHOLD){
			$bighash{"$gene1\t$gene2"} = $cor;
		}
	}	
	close IFIL;
}

print "\n\tREADING COMPLETE!\n\n\tCalculating MR scores now...\n";

open (OFIL, '>', $OUTPATH) or die "Couldn't write to $OUTPATH: $!\n";

my %rank;
my @decay_array = (5,10,25,50,100);
foreach my $pair (sort { $bighash{$b} <=> $bighash{$a} } keys %bighash) {
	#print "$pair\n";
	(my $gene1, my $gene2) = split(/\t/, $pair);
	my $cor = $bighash{$pair};
	
	$rank{$gene1}++;
	$rank{$gene2}++;
	
	my $mr = sqrt($rank{$gene1} * $rank{$gene2});
	print OFIL "$gene1\t$gene2\t$cor\t$mr";
	foreach my $decay (@decay_array){
		my $weight = exp(-1 * ($mr -1) / $decay);
		if ($weight < $WEIGHT_THRESHOLD){
			$weight = 0;
		}
		print OFIL "\t$weight";
	}
	print OFIL "\n";

}
close OFIL;


