#! /usr/bin/perl -w
open INFILE, "prefix.txt";
@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;

#create map file
#specific only for V4 primers, Stoeck et al.

for $i(@prefix)
{
   $filename = $i.'_map.txt';
   open OUTFILE, ">$filename";
   print OUTFILE '#',"SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n";
   $i =~ s/_//g;
   print OUTFILE $i, "\tCCAG\tCASCYGCGGTAATTCC\t", $i, "\n";
   close OUTFILE;
}
