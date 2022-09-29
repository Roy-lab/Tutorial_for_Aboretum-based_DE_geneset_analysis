#!/usr/bin/env perl
# Junha Shin, Sep 2022
# USAGE: ./transpose_matrix.pl [matrix.txt]
# Transpose a text matrix (tab-delimited) file WITH headers

use strict;
use warnings;
unless (scalar @ARGV == 1) {die " - USAGE: ./transpose_matrix.pl [matrix.txt]\n\n";}

my @rows = ();
my @cols = ();
my %data = ();
my $line_start = 0;
open (MAT, $ARGV[0]) or die " - Cannot open matrix file\n\n";
while (<MAT>) {
	my $line = $_; chomp $line;
	my @line = split /\t/, $line;

	if ($line_start==0) {
		shift @line;
		@rows = @line;
		$line_start++;
	} else {
		my $rowid = shift @line;
		push @cols, $rowid;
		for (my$i=0; $i<scalar @line; $i++) {
			$data{$rowid}{$rows[$i]} = $line[$i];
		}
	}
}
close MAT;

my @legend = @cols;
unshift @legend, "";
print join ("\t", @legend)."\n";
foreach my$row (@rows) {
	my @newline = ($row);
	foreach my$col (@cols) {
		push @newline, $data{$col}{$row};
	}
	print join ("\t", @newline)."\n";
}

exit		
