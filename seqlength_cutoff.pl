#!/usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;


my $file = $ARGV[0]; 
my $min = $ARGV[1];
my $max = $ARGV[2];
my $out = $ARGV[3];

open (FILE, ">>$out") or die ("Error : Cannot open file $out for writing..!\n");

my $seq_in  = Bio::SeqIO->new( -format => 'fasta',-file => $file);

while( my $seq1 = $seq_in->next_seq() ) {	
	
	my $id  = $seq1->primary_id;
	chomp $id;
	my $seq = $seq1->seq;
	chomp $seq;
	my $lseq = length($seq);
	if($lseq>=$min && $lseq <=$max){
		print FILE ">",$id,"\n",$seq,"\n";	
	}
}
