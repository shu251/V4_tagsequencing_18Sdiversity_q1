#! /usr/bin/perl -w
open INFILE, "prefix.txt";
@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;


#Perl script to automate V4 tag sequence QC process
#last updated 06-13-2017, Sarah K. Hu



#(1) Import prefix file, named "prefix.txt" and starting loop
for $i(@prefix)
{
#(2) Merge paired end (R1 and R2) reads with a 20 bp overlap
print "join_paired_ends.py -f ",$i,"_L001_R1_001.fastq -r ",$i,"_L001_R2_001.fastq -o joined_",$i," -j 20\n"; 
print "mv joined_",$i,"/fastqjoin.join.fastq ",$i,"_merged.fastq\n";

#(3) Perform split library command in QIIME, filter sequences with Q score 30 (-q 29)
print "split_libraries_fastq.py -i ",$i,"_merged.fastq -m ",$i,"_map.txt --barcode_type 'not-barcoded' --sample_ids ",$i," -q 29 -n 0 -o split_",$i,"\n";
print "mv split_",$i,"/seqs.fna ",$i,".merged.Q30.fasta\n";

#(4) Clip primers. Allows sequences to have either forward or reverse primers 
#(4.1) Searches for forward primer, discard seqs without those primers
print "cutadapt -g CCAGCASCYGCGGTAATTCC -O 3 --discard-untrimmed -m 10 -o ",$i,".assembled.clipped.regF.fasta ",$i,".merged.Q30.fasta >> ",$i,".filter.log\n";

#Name discarded reads from 4.1 "discard_regF.fasta"
print "cutadapt -g CCAGCASCYGCGGTAATTCC -O 3 --untrimmed-output discard_regF.fasta -m 10 -o tmp.fasta ",$i,".merged.Q30.fasta >> ",$i,".filter.log\n";
print "rm tmp.fasta\n"; #remove tmp (which is a duplicate)
#4.2 Check discarded for reverse primer
print "cutadapt -a TYRATCAAGAACGAAAGT -O 3 --discard-untrimmed -m 10 -o ",$i,".assembled.clipped.regR.fasta discard_regF.fasta>> ",$i,".filter.log\n";

#(5) Combine clipped sequences
print "cat ",$i,".assembled.clipped.regF.fasta ",$i,".assembled.clipped.regR.fasta >> ",$i,".assembled.clipped.fasta\n"; 

#(6) Move excess files to "split_*" directories which were created during the split library step
print "mv ",$i,"*reg* split_",$i,"/\n"; #move excess files from trimming to split file
print "mv ",$i,".filter.log split_",$i,"/\n";

# (7) Length filter: seqlength_cutoff.pl [input.fasta] [min] [max] [output.fasta] 
print "./seqlength_cutoff.pl ",$i,".assembled.clipped.fasta 150 500 ",$i,".assembled.clipped.len.fasta\n";

#(8) Chimera check with vsearch (uchime) using a reference database
print "vsearch --uchime_ref ",$i,".assembled.clipped.len.fasta --db /galadriel/sarah/PR2/pr2.qiime.fasta --uchimeout ",$i,".uchimeinfo_ref --chimeras ",$i,".chimeras_ref.fasta --strand plus --nonchimeras ",$i,".assembled.clipped.len.nc.final.fasta \n";

#(9) Move excess files (chimeras) to split directories
print "mv ",$i,".chimeras_ref.fasta split_",$i,"/\n"; #move excess chimera files to split dir
print "mv ",$i,".uchimeinfo_ref split_",$i,"/\n";

#(10) Combine all reads together
print "cat ",$i,".assembled.clipped.len.nc.final.fasta >> allseqs_test.fasta\n";

} 
