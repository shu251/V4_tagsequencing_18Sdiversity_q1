#! /usr/bin/perl -w
open INFILE, "prefix.txt";
@prefix = ();
while (<INFILE>)
{ chomp;
  push(@prefix, $_);
}
close INFILE;


#Perl script to automate V4 tag sequence QC process
#last updated 04-28-2018, Sarah K. Hu


for $i(@prefix)
{

#(2) Merge paired end (R1 and R2) reads with a 20 bp overlap
print "join_paired_ends.py -f raw_dir/",$i,"_L001_R1_001.fastq -r raw_dir/",$i,"_L001_R2_001.fastq -o excess_",$i," -j 20\n"; 
print "mv excess_",$i,"/fastqjoin.join.fastq ",$i,"_merged.fastq\n";

#(3) Perform split library command in QIIME, filter sequences with Q score 30 (-q 29)
print "split_libraries_fastq.py -i ",$i,"_merged.fastq -m map_dir/",$i,"_map.txt --barcode_type 'not-barcoded' --sample_ids ",$i," -q 29 -n 0 -o excess_",$i,"/split_",$i,"\n";
print "mv excess_",$i,"/split_",$i,"/seqs.fna ",$i,"_merged_QC.fasta\n";

#(4) Remove V4 primers
print "cutadapt -g CCAGCASCYGCGGTAATTCC ",$i,"_merged_QC.fasta > tmpFOR",$i,".fasta 2> excess_",$i,"/primerreportFOR.txt\n";

print "cutadapt -a TYRATCAAGAACGAAAGT tmpFOR",$i,".fasta > ",$i,"_merged_QC_trim.fasta 2> excess_",$i,"/primerreportREV.txt\n";

print "rm tmpFOR",$i,".fasta\n";

#(5) Length filter: seqlength_cutoff.pl [input.fasta] [min] [max] [output.fasta] 
print "./seqlength_cutoff.pl ",$i,"_merged_QC_trim.fasta 150 500 ",$i,"_merged_QC_trim_len.fasta\n";

#(6) Chimera check with vsearch (uchime de novo)
print "vsearch --uchime_denovo ",$i,"_merged_QC_trim_len.fasta --nonchimeras ",$i,"_merged_QC_trim_len_nc.fasta\n";

#(7) Get stats and clean up directory:
print "count_seqs.py -i raw_dir/",$i,"_L001_R1_001.fastq,",$i,"_merged.fastq,",$i,"_merged_QC.fasta,",$i,"_merged_QC_trim.fasta,",$i,"_merged_QC_trim_len.fasta,",$i,"_merged_QC_trim_len_nc.fasta > stats_",$i,".txt\n";

print "mv ",$i,"_merged.fastq excess_",$i,"\n";
print "mv ",$i,"_merged_QC.fasta excess_",$i,"\n";
print "mv ",$i,"_merged_QC_trim.fasta excess_",$i,"\n";
print "mv ",$i,"_merged_QC_trim_len.fasta excess_",$i,"\n";

#(8) Combine all reads together
print "cat ",$i,"_merged_QC_trim_len_nc.fasta >> allseqs_test.fasta\n";

}


