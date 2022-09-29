#!/usr/bin/env perl
#usage: ./script.pl
use strict;
use warnings;


my %mat = ();
open (F, "../clustermeans.txt") or die " - NO clustermeans.txt\n\n";
while (<F>) {
	my $line = $_; chomp $line;
	my @line = split /\t/, $line;
	my @temp = split /_/, $line[0];	#[0]=cols, [1]=clusterID
	$mat{$temp[1]}{$temp[0]} = $line[1];
}
close F;

my @order = ();
open (ORD, "../genemembers_perog.txt") or die;
#open (ORD, $ARGV[0]) or die;
while (<ORD>) {
	my $line = $_; chomp $line;
	my @line = split /\t/, $line;
	my @temp = split /;/, $line[1];
	foreach my$col (@temp) {
		unless ($col =~ /^Anc/) {
			my @t = split /_/, $col;
			my $colname = $t[0];
			push @order, $colname;
		}
	}
	last;
}
close ORD;

my $output = "order.txt";
open (O, "> $output") or die;
print O join ("\n", @order)."\n";
close O;


my @sorted = sort {$a <=> $b} keys %mat;

foreach my$id (@sorted) {
	my @newline = ();
	foreach my$om (@order) {
		push @newline, $mat{$id}{$om};
	}
	print join ("\t", @newline)."\n";
}

exit
