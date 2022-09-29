#!/usr/bin/env perl
# Junha Shin, Sep 2022
# USAGE: ./generating_meanvals.pl [matrix.txt] [output_prefix]
# INPUT matrix.txt		tab-delimited text file, [genes x cells] data matrix WITH headers
# INPUT output_prefix		string for naming gene IDs and output file, e.g. "c1"
# OUTPUT [prefix]_meanval.txt	FORMAT: col1=[prefix]_geneID, col2=mean values of all cell values


use strict;
use warnings;
unless (scalar @ARGV == 2) {die " - USAGE: ./generating_meanvals.pl [matrix.txt] [output_prefix]\n\n";}

my $input_file = $ARGV[0];
my $output_prefix = $ARGV[1]; chomp $output_prefix;
my $output_filename = $output_prefix."_meanval.txt";

my @gnames = ();
my $line_num = 0;
open (OUT, "> $output_filename") or die " - Cannot write output file.\n\n";
print " - Printing $output_filename..\n";
open (MAT, $input_file) or die " - Cannot open input matrix file.\n\n";
while (<MAT>) {
	my $line = $_; chomp $line;
	if ($line_num==0) {$line_num++; next;}
	my @line = split /\t/, $line;
	my $gnid = shift @line;
	$gnid = $output_prefix."_".$gnid;
	push @gnames, $gnid;

	my $sum = 0;
	foreach my$val (@line) {$sum += $val;}
	my $mean = $sum / scalar @line;
	$mean = sprintf "%.6f", $mean;
	
	print OUT $gnid."\t".$mean."\n";
}
close MAT;
close OUT;

exit
