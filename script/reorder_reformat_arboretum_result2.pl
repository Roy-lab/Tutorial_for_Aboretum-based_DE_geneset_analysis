#!/usr/bin/env perl
# Junha Shin, Sep 2022, revised Dec 2022
# USAGE: ./reorder_reformat_arboretum_result2.pl [arboretum_result_dir] [order.txt]

use strict;
use warnings;
unless (scalar @ARGV == 2) {die " - USAGE: ./reorder_reformat_arboretum_result2.pl [arboretum_result_dir] [order.txt]\n\n";}

my $resdir = $ARGV[0]; chomp $resdir;
open (ORD, $ARGV[1]) or die;
my @colorder = <ORD>; chomp @colorder;
close ORD;

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
my %sorted_position = ();
my $idx = 1;
foreach my$old_clid (@sorted) {
	$sorted_position{$old_clid} = $idx;
	$idx++;
}
## making old -> new ID map
my %change = ();	# k=old_clid, v=new_clid
my $new_clid = 0;
print "Before\tAfter\n";
foreach my$old_clid (@sorted) {
	$change{$old_clid} = $new_clid;
	print "$old_clid\t$new_clid\n";
	$new_clid++;
}
## writing out neworder.txt file
my $neworder = $resdir."/neworder.txt";
open (OUT, "> $neworder") or die;
my @before = sort {$a <=> $b} @sorted;
unshift @before, "before";
print join ("\t",@before)."\n";
print OUT join ("\t",@before)."\n";
my @newid = ("newID");
my %position_oldid = ();
foreach my$old_clid (@before) {
        if ($old_clid eq "before") {next;}
        push @newid, $sorted_position{$old_clid};
        $position_oldid{$sorted_position{$old_clid}} = $old_clid+1;
}
my @temp = sort {$a <=> $b} keys %position_oldid;
my @newordidx = ("newordidx");
foreach my$pidx (@temp) {
        push @newordidx, $position_oldid{$pidx};
}
print join ("\t", @newordidx)."\n";
print OUT join ("\t",@newordidx)."\n";
print join ("\t",@newid)."\n";
print OUT join ("\t",@newid)."\n";
close OUT;

## parsing clusterassign_multspecies.txt
my %data = ();
open (ASGN, "$resdir/clusterassign_multspecies.txt") or die " - Cannot open clusterassign_multspecies.txt\n\n";
while (<ASGN>) {
	my $line = $_; chomp $line;
	my @line = split /\t/, $line;
	my $gnid = shift @line;
	for (my$i=0; $i<scalar @line; $i++) {
		my @temp = split "=", $line[$i];
		my $colname = $temp[0];
		my $clid = $temp[1];
		$data{$gnid}{$colname} = $clid;
	#	print "$colname\t$clid\n";
	}
}
close ASGN;

## rewriting result as format for the findTransitionGenesets
my $output = $resdir."/allcelltypes_clusterassign_brk.txt";
open (OUT, "> $output") or die " - Cannot create output file\n\n";
my @legend = @colorder;
unshift @legend, "Loci";
print OUT join ("\t", @legend)."\n";
foreach my$gnid (sort keys %data) {
	my @newline = ($gnid);
	my $pass = 0;
	foreach my$coln (@colorder) {
	#	print "$gnid\t$coln\t$data{$gnid}{$coln}\n";
		if ($data{$gnid}{$coln} < 0) {$pass++;}
		push @newline, $change{$data{$gnid}{$coln}};
	}
	if ($pass > 0) {next;}
	print OUT join ("\t", @newline)."\n";
}
close OUT;


exit

