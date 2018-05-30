#### Last updated April 28, 2018
As with any bioinformatic pipeline, with each new dataset I analyze, this protocol evolves. While I plan to continue updating this repository, feel free to contact me with questions or comments - sarah.hu[at]usc.edu

## 18S rRNA gene tag-sequencing - V4
Below protocol describes approach used for 18S tag-sequencing, specifically from demultiplex paired end reads targetting the V4 hypervariable region. Below steps go from quality checking and filtering raw fastq reads to OTU clustering. This protocol is specific to using QIIME1 and employs subsampled open-reference OTU picking in QIIME. A full description of this can be found [here](http://qiime.org/tutorials/open_reference_illumina_processing.html).

### Required programs/files
* Contents of this repo, especially scripts: create.map.pl & seqlength_cutoff.pl
* [QIIME](http://qiime.org/install/install.html) -  v.1.9.1 or higher. Note as of Jan 2018 QIIME1 will not be supported. Future iterations of this protocol will use QIIME2.
* fastqjoin - this should install with QIIME. Alternatively you can use [PEAR](https://sco.h-its.org/exelixis/web/software/pear/doc.html#installing) to merge PE sequences.
* [cutadapt](http://cutadapt.readthedocs.io/en/stable/installation.html) v.1.15 or higher
* [vsearch](https://github.com/torognes/vsearch) - v1.11.1
* Reference database for downstream OTU clustering & taxonomy assignment. For 18S (microbial eukaryotic) work, I prefer [PR2](https://github.com/vaulot/pr2_database/wiki)
* R - see optional steps below for using R to generate stats on QC steps.

## Step 1 - Get started
To follow step by step instructions below, follow along with test files provided from here: zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1236641.svg)](https://doi.org/10.5281/zenodo.1236641)

Zenodo contents:
* test_fastq_small.zip - set of 4 test samples (R1 and R2) with 800 sequences each. Commands in steps below use Test01 from this dataset.
* test_fastq_large.zip - full set of raw sequences (R1 and R2 per sample) to provide the opportunity to run the below tutorial with a full dataset.

Otherwise, ensure raw fastq files are demultiplexed and not gzipped.

First, set up local directory with test files and generate a 'prefix.txt' file (used to create map files and automate process). Following steps only refer to "Test01" sample.
```
# Sort raw sequences (using the small dataset here)
unzip test_fastq_small.zip
mkdir raw_dir
mv Test*.fastq raw_dir/
# Generate map files (prefix.txt file required)
./create.map.pl
mkdir map_dir
mv *map.txt map_dir/
```
More information on map files: Map files in QIIME serve to de-multiplex reads (based on barcodes and indices) and are used in 'split_library' QIIME step. Here, I've split up the V4 forward primer to fulfill the necessary barcode / index space required in the QIIME mapping file. Since we get de-multiplexed fastq files from the sequence center directly, we do not need to input our barcode sequences here.
[Check mapping file](http://qiime.org/1.6.0/scripts/check_id_map.html)

 
## Step 2 - Merge paired end (R1 and R2) reads with a 20 bp overlap
Join the ~250 bp R1 and R2 to form an ~400 bp merged read. Below requires a 20 bp overlap.

```
join_paired_ends.py -f raw_dir/Test01_L001_R1_001.fastq -r raw_dir/Test01_L001_R2_001.fastq -o excess_Test01 -j 20
mv excess_Test01/fastqjoin.join.fastq Test01_merged.fastq
```

## Step 3 - Initial Q score filter
Perform split library command in QIIME, filter sequences with Q score 30
Here, the map file writes in the "Test01" name into the fasta header. Output of split_libraries_fastq.py generates a 'split' file with "seqs.fna", move and rename this file with 'mv'.

```
# Run quality check
split_libraries_fastq.py -i Test01_merged.fastq -m map_dir/Test01_map.txt --barcode_type 'not-barcoded' --sample_ids Test01 -q 29 -n 0 -o excess_Test01/split_Test01
# Move and rename
mv excess_Test01/split_Test01/seqs.fna Test01_merged_QC.fasta
```

## Step 4 - Primer removal
Commands below are specific for using the [V4 Stoeck et al. 2010 primers](https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1365-294X.2009.04480.x).
Searches for 5' primer, trims it and then searches for 3' primer. If a sequences does not have either primer, it is removed.
Generates a report file that summarizes rate of primer mismatches and number of reads removed.
```
# Remove forward primers, report in excess dir, generate a tmp file
# I'm allowing 30% error
cutadapt -g CCAGCASCYGCGGTAATTCC -e 0.30 --discard-untrimmed Test01_merged_QC.fasta > tmpFOR_Test01.fasta 2> excess_Test01/primerreportFOR_Test01.txt
# remove reverse primers, continue report, generate primer-free fasta file
cutadapt -a ACTTTCGTTCTTGATYRA -e 0.30 --discard-untrimmed tmpFOR_Test01.fasta > Test01_merged_QC_trim.fasta 2> excess_Test01/primerreportREV_Test01.txt
```

## Step 5 - Sequence length QC
Here remove sequences shorter than 150 bp and longer than 500 bps, and rename.
Usage: 
seqlength_cutoff.pl input.fasta min max output.fasta

```
./seqlength_cutoff.pl Test01_merged_QC_trim.fasta 150 500 Test01_merged_QC_trim_len.fasta
```

## Step 6 - Get QC statistics for each sample
This simply counts the number of sequences in each of the QCed files. Compiled stats*.txt file can be reviewed to see what steps removed the most (or fewest) total number of sequences. Depending on your research goals, it may be in your best interest to keep as much data as possible, i.e. relax the stringency during sequence qc.

```
count_seqs.py -i raw_dir/Test01_L001_R1_001.fastq,Test01_merged.fastq,Test01_merged_QC.fasta,Test01_merged_QC_trim.fasta,Test01_merged_QC_trim_len.fasta > stats_Test01.txt
# Clean up directory
mv Test01_merged.fastq excess_Test01
mv Test01_merged_QC.fasta excess_Test01
mv Test01_merged_QC_trim.fasta excess_Test01
mv Test01_merged_QC_trim_len.fasta excess_Test01
```

Output text files show you the total number of sequences at each step. For example, for "Test01", we started with 800 sequences, initial merging remove almost 200 sequences. Another 100 sequences were removed during the QC step. We've kept 69% of the sequences.
Now start an R environment to compile statistics of the QC run for each sample:

## Step 7 - Repeat for all sequences and combine all QCed reads together

To repeat, see automation pipeline information at [protocols.io](https://www.protocols.io/view/microbial-eukaryotic-18s-tag-sequence-processing-q-j9bcr2n) and using 'V4_tagseq_automate_QC_04282018.pl'

```
# Make sure to uncomment lines for chosen chimera checking step and replace "$PWD" with correct path.
perl V4_tagseq_automate_QC_04282018.pl > automate.sh
bash automate.sh
```

QIIME map file (during the split library step) introduced sample names (e.g. "Test01") into the fasta file header. So when we combine all of our fasta files together we can still track which sequences belong to which samples.

```
cat excess_*/Test*_trim_merged_Q30_len.fasta >> allseqs_test.fasta
```

## Step 8 - Chimera check and removal

I recommend pooled chimera checking, where chimera are removed from the entire dataset, rather than for each individual sample. There are two options for removing chimeric reads from the dataset, de novo or reference based. Both options are listed below. We're using [vsearch](https://github.com/torognes/vsearch) here.
```
# de novo chimera checking
vsearch --uchime_denovo allseqs_test.fasta --nonchimeras allseqs_test_denovoNC.fasta

# reference-based chimera checking
vsearch --uchime_ref allseqs_test.fasta --nonchimeras allseqs_test_refNC.fasta --db $PWD/[reference database] # insert path to reference database, e.g. PR2
```

## Step 9 - (optional) Get run statistics
After run is repeated (9a), use the below R script to import all of the .txt files generated in step 6 so we can look at how many sequences were removed at each QC step.
Since chimera checking (step 8) was done for the whole dataset, you'll have to add how many sequence were removed via chimera checking later.

Start R environment
```
# R
library(reshape2)
#
# Select files with "stats" in the name and end with ".txt"
stats_files <- intersect(list.files(pattern="stats"), list.files(pattern=".txt"))
# Run loop to compile all stats for each sample
for (file in stats_files){
  if (!exists("qc_info")){
    # Import w.o header. Ensure read separation is a space " "
    imported<-read.delim(file, header=FALSE, sep=" ")
    # Select columns and rows desired. Take note of order of QC steps:
    tmp<-imported[c(1,4)]
    qc_info<-tmp[1:5,]
    # Obtain sample ID and generate new column with sample ID
    split<- colsplit(qc_info$V4, "_", c("identity", "else"))
    samplename<-split[1,1]
    qc_info$SampleID<-samplename
    # Add in QC steps, *NOTE* will need to change order depending on how count_seqs.py was done
    qc_info$QC_step<-c("Primer removal", "Length trim", "Q score", "Merge", "Raw")
    colnames(qc_info)[1:2]<-c("Sequences", "File name")
    }
  if (exists("qc_info")){
    imported<-read.delim(file, header=FALSE, sep=" ")
    # Select columns and rows desired. Take note of order of QC steps:
    tmp<-imported[c(1,4)]
    qc_tmp<-tmp[1:5,] #instead - call qc_tmp
    # Obtain sample ID and generate new column with sample ID
    split<- colsplit(qc_tmp$V4, "_", c("identity", "else"))
    samplename<-split[1,1]
    qc_tmp$SampleID<-samplename
    # Add in QC steps, *NOTE* will need to change order depending on how count_seqs.py was done
    qc_tmp$QC_step<-c("Primer removal", "Length trim", "Q score", "Merge", "Raw")
    colnames(qc_tmp)[1:2]<-c("Sequences", "File name")
    qc_info<-rbind(qc_info, qc_tmp) # Now add to growing df called "qc_info"
    rm(qc_tmp) # Remove the tmp df created here
    qc_info<-unique(qc_info)
  }
}
head(qc_info)
# Cast to show by sample change
qc_stats<-dcast(qc_info[c(3:4,1)], SampleID~QC_step, value.var = "Sequences")
# Write to .csv file
write.csv(qc_stats, file="sequence_qc.csv")
```

## Step 10 - OTU clustering

The following are directions for running subsampled open reference OTU clustering in QIIME1 as described by [Rideout et al. 2014 PeerJ](https://peerj.com/articles/545/). 

Below steps use the output from the QC steps 1-10. 'allseqs_test.fasta' is also included in the repo so one can start here.

# Considerations for OTU clustering

* Depending on the questions you want to ask your data, decide about what type of OTU clustering you want to do (both algorithm and percent similarity).
* QIIME is a great resource for tutorials on [OTU clustering](http://qiime.org/tutorials/otu_picking.html#running-the-otu-picking-workflows)
* Obtain a reference database for these steps. For 18S (microbial eukaryotic) work, I prefer [PR2](https://github.com/vaulot/pr2_database/wiki)


## Step 10a - run clustering:

This command takes in the combined fasta file from above and outputs a directory with all necessary files. Replace "$PWD/pr2_DB.fasta" with the location of your reference database. Make sure you have an associated taxonomy file with the database. [pick_open_reference_otus.pu](http://qiime.org/scripts/pick_open_reference_otus.html)
There are several actions wrapped into this one command:
1. Closed reference OTU clustering is performed against the specified reference database. This step essentially separates your input sequences into those that did cluster ('step1_otus/_clusters.uc, _otus.txt, step1_rep_set.fna') and those that did NOT cluster ('step1_otus/_failures.txt, failures.fasta')
2. 'failures.fasta' is then randomly subsampled to create 'subsampled_failures.fna'. Percent subsample is specified based on '-s', default is 0.001 (0.1%).
3. Next, these 'failed' sequences are clustered de novo. Each cluster is labeled as "new reference sequence". This step generates 'step2_otus/subsampled_seqs_clusters.uc, subsampled_seqs_otus.txt, step2_rep_set.fna'
4. These newly generated de novo OTUs will now serve as the reference OTUs to perform closed reference OTU picking with the full 'failures.fasta' file (from step 1). 
5. Now you have the option to take those reads that were not clustered into OTUs or hit the reference database and conduct additional de novo OTU clustering. *The command written below suppresses this step ('--suppress_step4'). At this point, especially for microbial eukaryotes, I feel we have tried to 'salvage' as many of the reads as possible to detect potentially novel or underrepresented taxa. The final de novo OTU clustering step is more likely to produce OTUs that have very low sequence counts and will have no taxonomic assignments - thus we prefere to skip this step. '--suppress_taxonomy_assignment' is also passed, as we prefer to perform this step separately (allows for additional troubleshooting and testing). *NOTE* suppressing step 4 and taxonomic assignment will produce an error.
6. The output file 'final_otu_map_mc2.txt' is an OTU map comprised of OTUs that contained at least 2 sequences (global singletons removed!). Based on this OTU map, the 'rep_set.fna' file was generated.
```
#OTU clustering with uclust, will make a new directory (pick_open)
pick_open_reference_otus.py -i allseqs_test_refNC.fasta -o pick_open -m uclust -r $PWD/pr2_DB.fasta --suppress_step4 --suppress_taxonomy_assignment
```

## Step 10b - assign taxonomy using uclust
This step takes the representative sequences from the pick open reference clustering (10a) and assigns taxonomy using the PR2 (or db of choice). Reads are assigned taxonomic identities using uclust (-m). Successful assignment requires the hit to have at least 90% similarity (--similarity 0.90) and uclust will consider up to 3 hits (--uclust_max_accepts 3). [assign_taxonomy.py](http://qiime.org/scripts/assign_taxonomy.html)
Replace "$PWD/pr2_DB.fasta" and "$PWD/pr2_tax.fasta" with location of database and accompanying taxonomy file.
```
cd pick_open
# Assign taxonomy using PR2 database, creates a new directory (uclust_taxonomy)
assign_taxonomy.py -i rep_set.fna -t $PWD/pr2_tax.fasta -r $PWD/pr2_DB.fasta -m uclust --similarity 0.90 --uclust_max_accepts 3 -o uclust_taxonomy
```

## Step 10c - Generate OTU table
Generate an OTU table in a biom format. Then convert to a text file and include taxonomy information:

```
make_otu_table.py -i final_otu_map.txt -o V4_tagseq_test.biom -t uclust_taxonomy/rep_set_tax_assignments.txt 
biom convert -i V4_tagseq_test.biom -o V4_tagseq_test.txt --to-tsv --header-key=taxonomy
```

## Downstream analyses
See other repository: "PreliminaryFigures_V4_tagseq" for generating preliminary figures from OTU tables in R.
