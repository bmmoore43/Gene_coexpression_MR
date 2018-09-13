#!/usr/local/bin/perl
use strict;
use warnings;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);

###############################################
# Run: perl calc_clusterone_modules.pl -i [infile] -o [outfile] [options]
#
# calc_clusterone_modules.pl
#
#
# Jen Wisecaver
# 20180904, jwisecav@purdue.edu
################################################

# usage
my %opts; getopt('idcpq', \%opts );

my $usage = "\n\tUsage Error:
	Run: perl calc_clusterone_modules.pl -i [infile] -o [outfile] -c [clusterone] -d [decay] [options]
	
	OPTIONS =
	Files:
	-i: <FILENAME> INFILE: path to input file containing mutual ranks
	-c: <FILENAME> CLUSTERONE: path to clusterone jar file
	-d: <INTEGER> DECAY: decay rate to use for clusterone network (currently can be 5, 10, 25, 50, 100)
	
	General:
	-p: <NUMBER> PVAL_THRESHOLD: retained modules must have p value <= this value (default = 1 to retain all modules)
	-q: <NUMBER> QUAL_THRESHOLD: retained modules must have weight >= this value (default = 0 to retain all modules)

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

my $DECAY = 0; 
if ($opts{d}) { if ($opts{d} =~ m/^\d+$/) { 
	if (int($opts{d}) eq $opts{d}) { $DECAY = "$opts{d}"; }} else { die "\n\tInput <d>\n$usage"; } 
}else{ die "\n\tInput <d>\n$usage"; }

my $CLUSTERONE;
if ($opts{c}) { 
	$CLUSTERONE = "$opts{c}";
}else{ die "\n\tInput <c>$usage"; }
my $program="java -jar $CLUSTERONE";

my $PVAL_THRESHOLD = 1; 
if ($opts{p}) { if ( looks_like_number($opts{p}) ) { $PVAL_THRESHOLD = "$opts{p}"; } else { die "\n\tInput <p>\n$usage"; } }

my $QUAL_THRESHOLD = 0; 
if ($opts{'q'}) { if ( looks_like_number($opts{'q'}) ) { $QUAL_THRESHOLD = "$opts{'q'}"; } else { die "\n\tInput <q>\n$usage"; } }


#############################
###### CREATE ABC FILE ######
#############################
open (IFIL, '<', $INPATH) or die "Couldn't open file $INPATH: $!\n";

(my $abcfile = $INFILE) =~ s/\.txt$//;
$abcfile = $abcfile . $DECAY . '.abc';
(my $csvfile = $abcfile) =~ s/\.abc$/.modules.csv/;
(my $outfile = $abcfile) =~ s/\.abc$/.modules.txt/;

my %decay_hash = (5, 4, 10, 5, 25, 6, 50, 7, 100, 8);

open (OFIL, '>', $abcfile) or die "Couldn't write to file $abcfile: $!\n";
while ( my $line = <IFIL> ) {
	chomp $line;

	my @col = split(/\t/, $line);
	my $gene1 = $col[0];
	my $gene2 = $col[1];
	my $weight = $col[$decay_hash{$DECAY}];
	next if ($weight == 0);
	
	print OFIL "$gene1\t$gene2\t$weight\n";
}
close IFIL;
close OFIL;

############################
###### RUN CLUSTERONE ######
############################
my $command = "$program $abcfile --output-format csv > $csvfile";
my $status = system($command); 
if ($status != 0) { print "clusterone failed! COMMAND: $command\n"; }

##############################
###### PARSE CLUSTERONE ######
##############################
open (IFIL, '<', $csvfile) or die "Couldn't open file $csvfile: $!\n";
open (OFIL, '>', $outfile) or die "Couldn't write to file $outfile: $!\n";
while ( my $line = <IFIL> ) {
	chomp $line;
	next if ($line =~ /^Cluster/);
	my @col = split(/,/, $line);

	my $curr_clust = $col[0];
	$curr_clust = sprintf("%06d", $curr_clust);
	$curr_clust = 'NET' . $DECAY . 'MOD' . $curr_clust;
	
	my $quality = $col[5];
	my $pvalue = $col[6];
	next if ($pvalue > $PVAL_THRESHOLD);
	next if ($quality < $QUAL_THRESHOLD);
	
	my $gene_list = $col[7];
	$gene_list =~ s/"//g;
	my @gene_array = split(/ /, $gene_list);
	
	print OFIL "$curr_clust\t$quality\t$pvalue\t@gene_array\n";
}
close IFIL;
close OFIL;
