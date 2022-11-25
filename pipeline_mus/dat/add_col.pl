#! /usr/bin/perl

use strict; use warnings;

open(FILE,'<',"mouse_aa.txt") or die "ERRaa\n"; my @aa=<FILE>; chomp @aa; close(FILE);
open(FILE,'<',"mouse_dat.txt") or die "ERRdt\n"; my @dt=<FILE>; chomp @dt; close(FILE);

my (%aa, %dt); my @prot;
for my $line (@aa) { my @a=split("\t",$line); $aa{$a[0]} = $line; push(@prot,$a[0]) }
for my $line (@dt) { my @a=split("\t",$line); $dt{$a[0]} = $line }  

for my $p (@prot) {
  if ( length($aa{$p}) && length($dt{$p}) ) {
    print "$dt{$p}";
    my @a = split("\t", $aa{$p});
    for my $ind (1 .. $#a) { print "\t$a[$ind]" }
    print "\n";
  }
  else { print "!! Missing $p\n" }
}
