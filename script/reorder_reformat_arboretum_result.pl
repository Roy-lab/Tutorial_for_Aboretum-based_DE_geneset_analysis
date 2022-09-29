#!/usr/bin/env perl
# Junha Shin, Sep 2022
# USAGE: ./reorder_reformat_arboretum_result.pl [arboretum_result_dir]

use strict;
use warnings;
unless (scalar @ARGV == 1) {die " - USAGE: ./reorder_reformat_arboretum_result.pl [arboretum_result_dir]\n\n";}

my $resdir = $ARGV[0]; chomp $resdir;

## taking sum of all cluster mean valuse
my %mat = ();
open (F, "$resdir/clustermeans.txt") or die " - Cannot open clustermeans.txt\n\n";
while (<F>) {
        my $line = $_; chomp $line;
        my @line = split /\t/, $line;
        my @temp = split /_/, $line[0]; #[0]=cols, [1]=clusterID
        #$mat{$temp[1]}{$temp[0]} = $line[1];
        $mat{$temp[1]} += $line[1];	# summation of cluster mean valuess of samples
}
close F;

## ordering by sum value
my @sorted = sort {$mat{$a} <=> $mat{$b}} keys %mat;	# ascending order

## making old -> new ID map
my %change = ();	# k=old_clid, v=new_clid
my $new_clid = 0;
print "Before\tAfter\n";
foreach my$old_clid (@sorted) {
	$change{$old_clid} = $new_clid;
	print "$old_clid\t$new_clid\n";
	$new_clid++;
}

## parsing column orders for the head line printing
my @order = ();
open (ORD, "$resdir/genemembers_perog.txt") or die " - Cannot open genemembers_perog.txt\n\n";
while (<ORD>) {
        my $line = $_; chomp $line;
        my @line = split /\t/, $line;
        my @temp = split /;/, $line[1];
        foreach my$col (@temp) {
                unless ($col =~ /^Anc/) {
                        my @t = split /_/, $col;
                        my $colname = $t[0];
                        push @order, $colname;
		} else {
			push @order, $col;
		}
        }
        last;
}
close ORD;

## rewriting result as format for the findTransitionGenesets
my $output = $resdir."/allcelltypes_clusterassign_brk.txt";
open (OUT, "> $output") or die " - Cannot create output file\n\n";
open (RES, "$resdir/allspecies_clusterassign_lca_brk.txt") or die " - Cannot open allspecies_clusterassign_lca_brk.txt\n\n";
my @legend = @order;
unshift @legend, "Loci";
print OUT join ("\t", @legend)."\n";
while (<RES>) {
	my $line = $_; chomp $line;
	if ($line =~ /^Dummy/) {
		print OUT $line."\n";
	} else {
		my @line = split /\t/, $line;
		my $gnID = shift @line;
		my @newline = ($gnID);
		foreach my$clid (@line) {
			push @newline, $change{$clid};
		}
		print OUT join ("\t", @newline)."\n";
	}
}
close RES;
close OUT;

exit

