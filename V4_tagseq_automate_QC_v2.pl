#! /usr/bin/perl -w
open INFILE, "prefix.txt";
@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;


#Perl script to automate V4 tag sequence QC process
#last updated 09-29-2017, Sarah K. Hu


for $i(@prefix)
{
#(1) Trim primers with trimmomatic
print "java -jar /usr/local/bioinf/Trimmomatic-0.32/trimmomatic-0.32.jar PE ",$i,"_L001_R1_001.fastq ",$i,"_L001_R2_001.fastq ",$i,"_trim_R1_PE.fastq ",$i,"_trim_R1_orphan.fastq ",$i,"_trim_R2_PE.fastq ",$i,"_trim_R2_orphan.fastq LEADING:10 TRAILING:10 SLIDINGWINDOW:10:30 MINLEN:50 ILLUMINACLIP:trimPE_V4primer.fasta:2:30:10\n";

#(2) Merge paired end (R1 and R2) reads with a 20 bp overlap
print "join_paired_ends.py -f ",$i,"_trim_R1_PE.fastq -r ",$i,"_trim_R2_PE.fastq -o out_",$i," -j 20\n"; 
print "mv out_",$i,"/fastqjoin.join.fastq ",$i,"_trim_merged.fastq\n";

#(3) Perform split library command in QIIME, filter sequences with Q score 30 (-q 29)
print "split_libraries_fastq.py -i ",$i,"_trim_merged.fastq -m ",$i,"_map.txt --barcode_type 'not-barcoded' --sample_ids ",$i," -q 29 -n 0 -o split_",$i,"\n";
print "mv split_",$i,"/seqs.fna ",$i,"_trim_merged_Q30.fasta\n";
print "mv split_",$i,"/ out_",$i,"/ \n";

#(4) Length filter: seqlength_cutoff.pl [input.fasta] [min] [max] [output.fasta] 
print "./seqlength_cutoff.pl ",$i,"_trim_merged_Q30.fasta 150 500 ",$i,"_trim_merged_Q30_len.fasta\n";

#(5) Chimera check with vsearch (uchime) using a reference database
print "vsearch --uchime_ref ",$i,"_trim_merged_Q30_len.fasta --db /galadriel/sarah/PR2/pr2.qiime.fasta --uchimeout ",$i,".uchimeinfo_ref --chimeras ",$i,".chimeras_ref.fasta --strand plus --nonchimeras ",$i,"_trim_merged_Q30_len_nc.fasta \n";

#(9) Move excess files (chimeras) to split directories
print "mv ",$i,".chimeras_ref.fasta out_",$i,"/ \n";
print "mv ",$i,".uchimeinfo_ref out_",$i,"/ \n";
print "mv ",$i,"*orphan.fastq out_",$i,"/ \n";

#(10) Combine all reads together
print "cat ",$i,"_trim_merged_Q30_len_nc.fasta >> allseqs_test.fasta\n";

}

