#### Last updated September 30, 2017
As with any bioinformatic pipeline, with each new dataset I analyze, I evolve this protocol. While I plan to continue updating this repository, feel free to contact me with questions or comments - sarah.hu[at]usc.edu

## V4 tag sequence QC automation
Protocol describes quality checking process for raw fastq sequences. Includes some discussion for alternate approaches. End product is ready for OTU clustering. I've also included suggestions for OTU clustering at the end (specific for 18S data). Much of this pipeline can be applied for any tag sequencing analysis, but keep in mind this was specifically constructed to analyze single-celled microbial eukaryotic communities.

### Getting started
Step by step instructions below. Included perl script here automates the entire process. Directions on automating are included at [protocols.io](https://www.protocols.io/view/microbial-eukaryotic-18s-tag-sequence-processing-q-j5acq2e). For first timers, we suggest going through all steps below with one set of fastq samples (e.g. R1 and R2s) to familiarize yourself with process. 

This protocol (and my preferred method) is to perform as much quality checking on each individual sample before combining for downstream OTU clustering (or other analysis). This way, I can track any samples that may be not ideal for downstream analysis.


### Prerequisites
Programs/files required:
* Contents of this repo
* Database (preferrably with associated taxa list) to align sequences to. A version of the PR2 database formatted for QIIME can be found [here](https://drive.google.com/drive/folders/0Bxzw_UrYS4IAaEQ0b0lMQ0ZiUkU?usp=sharing)
* [QIIME](http://qiime.org/install/install.html) -  v.1.9.1 or higher. Note as of Jan 2018 QIIME1 will not be supported. Future iterations of this protocol will use QIIME2.
* fastqjoin - this should install with QIIME. Alternatively you can use [PEAR](https://sco.h-its.org/exelixis/web/software/pear/doc.html#installing) to merge PE sequences.
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) - v.0.32
* [vsearch](https://github.com/torognes/vsearch) - v1.11.1
* abyss (optional) - see step 9
* Make sure you know your primer sequences and generate a fasta file: 'trimPE_V4primer.fasta'. In the below example I'm using [V4 Stoeck et al. 2010 primers](http://onlinelibrary.wiley.com.libproxy1.usc.edu/doi/10.1111/j.1365-294X.2009.04480.x/full/) 
* Create a "prefix.txt" file to index sample fastq files.
* Run 'create.map.pl'
```
#Required the prefix.txt file
./create.map.pl
```
Map files: Map files in QIIME serve to de-multiplex reads (based on barcodes and indices) and are used in 'split_library' QIIME step. Here, I've split up the V4 forward primer to fulfill the necessary barcode / index space required in the QIIME mapping file. Since we get de-multiplexed fastq files from the sequence center directly, we do not need to input our barcode sequences here.
[Check mapping file](http://qiime.org/1.6.0/scripts/check_id_map.html)

### Step by step instructions

Following steps only refer to "Test01" sample.

## Step 1 - Remove primers using Trimmomatic

Removes primers. Imports a file that lists your forward and reverse sequences. See 'trimPE_V4primer.fasta'

* LEADING:10 TRAILING:10 - each base requires minimum quality of Q10 in leading and trailing ends
* SLIDINGWINDOW:10:30 - In each 10bp sliding window, the average quality required is Q30
* MINLEN:100 - All reads needs to be at least 50 bps in length
* Remove the V4 primers (listed in fasta file), allow 2 mismatches

```
java -jar /usr/local/bioinf/Trimmomatic-0.32/trimmomatic-0.32.jar PE Test01_L001_R1_001.fastq Test01_L001_R2_001.fastq Test01_trim_R1_PE.fastq Test01_trim_R1_orphan.fastq Test01_trim_R2_PE.fastq Test01_trim_R2_orphan.fastq LEADING:10 TRAILING:10 SLIDINGWINDOW:10:30 MINLEN:100 ILLUMINACLIP:trimPE_V4primer.fasta:2:30:10
```
 
## Step 2 - Merge paired end (R1 and R2) reads with a 20 bp overlap

```
join_paired_ends.py -f Test01_trim_R1_PE.fastq -r Test01_trim_R2_PE.fastq -o out_Test01 -j 20

mv out_Test01/fastqjoin.join.fastq Test01_trim_merged.fastq
```

## Step 3 - Initial Q score filter
Perform split library command in QIIME, filter sequences with Q score 30
Here, the map file writes in the "Test01" name into the fasta header. 

```
split_libraries_fastq.py -i Test01_trim_merged.fastq -m Test01_map.txt --barcode_type 'not-barcoded' --sample_ids Test01 -q 29 -n 0 -o split_Test01
mv split_Test01/seqs.fna Test01_trim_merged_Q30.fasta
mv split_Test01/ out_Test01/
```

## Step 4 - Length filter
Here remove sequences shorter than 150 bp and longer than 500 bps, and rename.
Usage: 
seqlength_cutoff.pl input.fasta min max output.fasta

```
./seqlength_cutoff.pl Test01_trim_merged_Q30.fasta 150 500 Test01_trim_merged_Q30_len.fasta
```

## Step 5 - Chimera check 
Remove chimeras using [vsearch](https://github.com/torognes/vsearch) (uchime) and the PR2 database
You will need to acquire a reference database that is best for your sample type. I am using the PR2 database. Alternatively, you can change the command to run vsearch chimera checking *de novo*. This will take longer.
```
vsearch --uchime_ref Test01_trim_merged_Q30_len.fasta --db /galadriel/db/PR2/pr2.qiime.fasta --uchimeout Test01.uchimeinfo_ref --chimeras Test01.chimeras_ref.fasta --strand plus --nonchimeras Test01_trim_merged_Q30_len_nc.fasta 

#Move excess files to out directories
mv Test01.chimeras_ref.fasta out_Test01/
mv Test01.uchimeinfo_ref out_Test01/
mv Test01*orphan.fastq out_Test01/
```

## Step 6 - Repeat for all sequences & get stats (optional)

To repeat, see automation pipeline information at [protocols.io](dx.doi.org/10.17504/protocols.io.g33byqn) and using 'V4_tagseq_automate_QC.pl'

To get sequence count information, you can use count_seqs.py in QIIME. I recommend running some kind of statistics on fasta files from each step. This way you can know exactly why a sample "failed" QC. Alternatively, you can see which steps in your pipeline are more strict (remove more sequences) or relaxed (didn't remove a lot of sequences). 

```
count_seqs.py -i Test01_L001_R1_001.fastq,Test01_L001_R2_001.fastq,Test01_trim_R1_PE.fastq,Test01_trim_R2_PE.fastq,Test01_trim_merged.fastq,Test01_trim_merged_Q30.fasta,Test01_trim_merged_Q30_len.fasta,Test01_trim_merged_Q30_len_nc.fasta
#Output:
315  : Test01_trim_merged_Q30_len_nc.fasta (Sequence lengths (mean +/- std): 419.0127 +/- 3.3219)
384  : Test01_trim_merged_Q30.fasta (Sequence lengths (mean +/- std): 418.8750 +/- 3.2612)
384  : Test01_trim_merged_Q30_len.fasta (Sequence lengths (mean +/- std): 418.8750 +/- 3.2612)
390  : Test01_trim_merged.fastq (Sequence lengths (mean +/- std): 415.9846 +/- 23.9856)
608  : Test01_trim_R2_PE.fastq (Sequence lengths (mean +/- std): 203.2648 +/- 34.6430)
608  : Test01_trim_R1_PE.fastq (Sequence lengths (mean +/- std): 231.1234 +/- 36.6052)
800  : Test01_L001_R2_001.fastq (Sequence lengths (mean +/- std): 249.9512 +/- 3.7032)
800  : Test01_L001_R1_001.fastq (Sequence lengths (mean +/- std): 250.6350 +/- 3.6748)
4289  : Total
```
This output shows you the total number of sequences at each step. We started with 800 sequences, initial quality trimming and primer removal/QC removed almost 200 sequences. Merging and additional QC removed 224 sequences, while length filtration removed no additional sequences. Finally, the chimera checking step only removed 69 sequences. We only kept around 40% of the total sequences. See next steps (below) for thoughts on when to remove samples.

If you have installed abyss or have another way to get quick stats on fasta files, use it here!
It is important to track how many sequences were lost during the QC process. Run this bit of code to obtain separate .txt files for each step.
```
abyss-fac -v *_L001_R1_001.fastq >>stats_raw.txt
abyss-fac -v *_trim_R1_PE.fastq >>stats_primerclipped.txt
abyss-fac -v *_trim_merged >>stats_merged.txt
abyss-fac -v *_trim_merged_Q30.fasta >> stats_Q30.txt
abyss-fac -v *_trim_merged_Q30_len.fasta >> stats_length_cutoff.txt
abyss-fac -v *_trim_merged_Q30_len_nc.fasta >> stats_chimeras.txt
```
Text file outputs from this can then be compiled using an R script to generate a QC stats file.
```
#Start R
stats<-c("stats_raw.txt", "stats_primerclipped.txt","stats_merged.txt", "stats_Q30.txt",  "stats_length_cutoff.txt", "stats_chimeras.txt")

for (file in stats){
  if (!exists("dataset")){
    dataset<-read.delim(file, header=TRUE, sep="\t",row.names=NULL)
    dataset<-dataset[c(2,4,8,10)]
    colnames(dataset)[1:4]<-c(file,"Min (bp)", "Max (bp)", "File name")
  }
  if (exists("dataset")){
    tmpdata<-read.delim(file, header=TRUE, sep="\t",row.names=NULL)
    tmpdata<-tmpdata[c(2,4,8,10)]
    colnames(tmpdata)[1:4]<-c(file,"Min (bp)", "Max (bp)", "File name")
    dataset<-cbind(dataset, tmpdata)
    rm(tmpdata)
  }
}
head(dataset)
write.csv(dataset, file="Seq_stats_QC.csv")
```

## Step 7 - Combine all reads together. QIIME map file (during the split library step) introduced sample names (e.g. "Test01") into the fasta file header. So when we combine all of our fasta files together we can still track which sequences belong to which samples.

```
cat Test*_trim_merged_Q30_len_nc.fasta >> allseqs_test.fasta
```

See related [protocols.io](https://www.protocols.io/view/microbial-eukaryotic-18s-tag-sequence-processing-q-j5acq2e) for directions on automating this process.

### Next steps - OTU clustering
I've included the test file result "allseqs_test.fasta" here as well.

# Considerations for OTU clustering
Depending on the questions you want to ask your data, decide about what type of OTU clustering you want to do (both algorithm and percent similarity).
 
QIIME is a great resource for tutorials on [OTU clustering](http://qiime.org/tutorials/otu_picking.html#running-the-otu-picking-workflows)

In this example I will use open reference OTU picking (Rideout et al. 2014, PeerJ), via QIIME and the PR2 database.
```
#OTU clustering with uclust, will make a new directory (pick_open)
pick_open_reference_otus.py -i allseqs_test.fasta -o pick_open -m uclust -r /galadriel/sarah/PR2/pr2.qiime.fasta --suppress_step4 --suppress_taxonomy_assignment
cd pick_open
#assign taxonomy using PR2 database, creates a new directory (uclust_taxonomy)
assign_taxonomy.py -i rep_set.fna -t /galadriel/db/PR2/ids.names.2.txt -r /galadriel/db/PR2/pr2.qiime.fasta -m uclust -o uclust_taxonomy
#make an OTU table, and convert to txt
make_otu_table.py -i final_otu_map.txt -o V4_tagseq_test.biom -t uclust_taxonomy/rep_set_tax_assignments.txt 
biom convert -i V4_tagseq_test.biom -o V4_tagseq_test.txt --to-tsv --header-key=taxonomy
```
### Use a text editor to remove the comment from the second line
Keep the "Constructed from biom file" commented out line.
See other repository: "PreliminaryFigures_V4_tagseq" for generating preliminary figures from OTU tables.

## Contributing
With helpful guidance and input from Jay Liu.
