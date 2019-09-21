#!/usr/bin/perl
use warnings;
use Getopt::Std;

###############################################
# Run: perl calc_pearsons_correlation.pl -i [infile] -o [outfile] [options]
#
# calc_pearsons_correlation.pl
#
#
# Jen Wisecaver
# 20180904, jwisecav@purdue.edu
################################################

# usage
my %opts; getopt('iot', \%opts );

my $usage = "\n\tUsage Error:
	Run: perl calc_pearsons_correlation.pl -i [infile] -o [outfile] [options]
	
	OPTIONS =
	Files:
	-i: <FILENAME> INFILE: path to input file: matrix of transformed & filtered gene counts (Header should start with a 'tab' see example_matrix.txt input file)
	-o: <DIRECTORY> OUTDIR: path to output directory
	
	General:
	-t: <INTEGER> THREADS: multithreaded using perl fork (default = 1)

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

my $OUTPATH;
if ($opts{o}) { 
	$OUTPATH = "$opts{o}";
}
else { die "\n\tInput <o>$usage"; }
unless (-d $OUTPATH) {
    print "making $OUTPATH";
    system("mkdir $OUTPATH");
}

my $THREADS = 1; 
if ($opts{t}) { if ($opts{t} =~ m/^\d+$/) { 
	if (int($opts{t}) eq $opts{t}) { $THREADS = "$opts{t}"; }} else { die "\n\tInput <t>\n$usage"; } }

##########################
###### Parse matrix ######
##########################
open (IFIL, '<', $INPATH) or die "Couldn't open file $INPATH: $!\n";
print "Reading in gene expression values from $INFILE...\n\n";
my %gene_hash;
my $pass = 0;
my $fail = 0;
my @array = <IFIL>;
shift @array;
foreach my $line (@array) {
	chomp $line;
	my @col = split(/\t/, $line);
	my $gene = shift @col;
	$gene =~ s/"//g;
	#print "$gene\n";
		
	@{$gene_hash{$gene}} = @col;
}
close IFIL;
print "\nCalculating Pearson's Correlations...\n\n";

my $nbdata = 0;
my @x;

my @childs; 
my @gene_array = keys %gene_hash;
my $num_genes = @gene_array;
my $subcount = $num_genes/$THREADS; #divide total number of genes by number of threads
$subcount=int(++$subcount); #round up to nearest whole number
for ( my $i = $THREADS; $i >= 1; $i--) {

    # extract a sub list of gene names from original file glob array
    my @subgenes = splice(@gene_array,0,$subcount);

    # convert to string because subroutines can only accept scalars.
    my $substring = join ("\t", @subgenes);

    my $pid = fork();
    if ($pid) {		# this is the parent
    	push(@childs, $pid);
    } elsif ($pid == 0) {	# this is a child
        &sub_fork(@subgenes);	# send substring to the subroutine 
        exit 0;		# this is important
    } else {
        die "couldnt fork: $!\n";
    } 
}

foreach (@childs) {
    my $tmp = waitpid($_, 0);
    print "\tdone with pid $tmp\n"; 
}

print "\nAll gene correlations printed to $OUTPATH\n";

sub sub_fork {
	my @sub_gene_array = @_;
	
# 	my $gene_count = @sub_gene_array;
# 	print "$gene_count in fork\n";
	for my $gene1 (@sub_gene_array){
		my $i1 = 0;
		my $sum1 = 0;
		for my $num (@{$gene_hash{$gene1}}){
			$i1++;
			$x[1][$i1]=$num;
			$sum1 = $sum1 + $num;
		}
		$nbdata=$i1;

		my %sort;
		for my $gene2 (keys %gene_hash){

			my $i2 = 0;
			my $sum2 = 0;
			for my $num (@{$gene_hash{$gene2}}){
				$i2++;
				$x[2][$i2]=$num;
				$sum2 = $sum2 + $num;
			}

			my $cor = &Correlation;
			#print "$gene1\t$gene2\t$cor\n";
			$sort{$gene2} = $cor;
		}

		my $outfile = "$OUTPATH/$gene1";
		open (OFIL, ">$outfile") or die "Cannot write to $outfile\n";
		foreach my $gene2 (sort { $sort{$b} <=> $sort{$a} } keys %sort){
			print OFIL "$gene2\t$sort{$gene2}\n";
		}
		close OFIL;
	}
}

sub Correlation {

$mean[1]=&Mean(1);
$mean[2]=&Mean(2);

$ssxx=&SS(1,1);
$ssyy=&SS(2,2);
$ssxy=&SS(1,2);

$correl=&Correl($ssxx,$ssyy,$ssxy);

$xcorrel=sprintf("%.4f",$correl);

return $xcorrel;

}  # End of Correlation


##########################################################################################
### Mean

sub Mean {

my ($a)=@_;
my ($i,$sum)=(0,0);

for ($i=1;$i<=$nbdata;$i++){
  $sum=$sum+$x[$a][$i];
}
$mu=$sum/$nbdata;

return $mu;

}

##########################################################################################
### SS = sum of squared deviations to the mean

sub SS {

my ($a,$b)=@_;
my ($i,$sum)=(0,0);

for ($i=1;$i<=$nbdata;$i++){
  $sum=$sum+($x[$a][$i]-$mean[$a])*($x[$b][$i]-$mean[$b]);
}

return $sum;

}

##########################################################################################
### Correlation

sub Correl {

my($ssxx,$ssyy,$ssxy)=@_;

$sign=$ssxy/abs($ssxy);

$correl=$sign*sqrt($ssxy*$ssxy/($ssxx*$ssyy));

return $correl;

}
