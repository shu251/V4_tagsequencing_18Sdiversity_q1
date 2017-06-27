#### As of June 2017
This project is under construction! It is a work in progress. If you would like more information, contact me directly at sarah.hu@usc.edu
Also, please feel free to contact me with suggestions or feedback on this pending repo

## V4 tag sequence QC automation
Protocol describes quality checking process for raw fastq sequences. Includes some discussion for alternate approaches. End product is ready for OTU clustering.

I've also included suggestions for OTU clustering at the end (specific for 18S data).

### Getting started
Step by step instructions below. Included perl script here automates the entire process. Directions on automating are included here (protocols.io) link. For first timers, we suggest going through all steps below with one set of fastq samples (e.g. R1 and R2s) to familiarize yourself with process. 

This protocol (and my preferred method) is to perform as much quality checking on each individual sample before combining for downstream OTU clustering (or other analysis). This way, I can track any samples that may be not ideal for downstream analysis.

To follow along with below step by step instructions:
- Download and unzip test_fastq_files.zip
- Download other scripts: create.map.pl, seqlength_cutoff.pl
- Download PR2 database from Caron Lab hosted drive: https://drive.google.com/drive/folders/0Bxzw_UrYS4IAaEQ0b0lMQ0ZiUkU?usp=sharing


### Prerequisites
Programs required:

QIIME - I'm using v.1.9.1
Install here: http://qiime.org/install/install.html
Make sure fastqjoin is installed. Alternatively you can use PEAR to merger PE sequences: https://sco.h-its.org/exelixis/web/software/pear/doc.html#installing

cutadapt - I'm using v.1.10: http://cutadapt.readthedocs.io/en/stable/guide.html

vsearch - v1.11.1: https://github.com/torognes/vsearch

abyss (optional) - see step 9

To know about your specific samples:
- Forward and reverse primer sequences, I'm using V4 Stoeck et al. 2010 primers 
- Preferred reference database, here, I'm using the PR2 reference database

## Automating each step

Before automating this protocol for your samples, go through the test data starting at Step 1 below. 

To automate run, first make sure the prefix.txt file has been made for your samples. See example.
Run the perl script V4_tagseq_automate_QC.pl
```
perl V4_tagseq_automate_QC.pl > run_qc.sh
```
Check the output and run!

### Step by step instructions
## Step 1 - Create map file
Map files in QIIME serve to de-multiplex reads (based on barcodes and indices) and are used in 'split_library' QIIME step. Here, I've split up the V4 forward primer to fulfill the necessary barcode / index space required in the QIIME mapping file. However, we aren't using this feature of the mapping file. 
 
Use this QIIME command to check mapping file:
http://qiime.org/1.6.0/scripts/check_id_map.html

```
./create.map.pl
```
This reads the prefix.txt file to create map files specific to the test fastq files for this example.

Following steps only refer to "Test01" sample.


## Step 2 - Merge paired end (R1 and R2) reads with a 20 bp overlap

```
join_paired_ends.py -f Test01_L001_R1_001.fastq -r Test01_L001_R2_001.fastq -o joined_Test01 -j 20

mv joined_Test01/fastqjoin.join.fastq Test01_merged.fastq
```

## Step 3 - Initial Q score filter
Perform split library command in QIIME, filter sequences with Q score 30

```
split_libraries_fastq.py -i Test01_merged.fastq -m Test01_map.txt --barcode_type 'not-barcoded' --sample_ids Test01 -q 29 -n 0 -o split_Test01
mv split_Test01/seqs.fna Test01.merged.Q30.fasta
```

## Step 4 - Remove primers

Clip V4 primers. Allows sequences to have either forward or reverse primers 
Note I am using Stoeck et al. 2010 V4 primers specifically for microbial eukaryotes

```
#look for forward primer
cutadapt -g CCAGCASCYGCGGTAATTCC -O 3 --discard-untrimmed -m 10 -o Test01.assembled.clipped.regF.fasta Test01.merged.Q30.fasta >> Test01.filter.log

#generate fasta file for seqs without a forward primer
cutadapt -g CCAGCASCYGCGGTAATTCC -O 3 --untrimmed-output discard_regF.fasta -m 10 -o tmp.fasta Test01.merged.Q30.fasta >> Test01.filter.log
rm tmp.fasta

#search those sequences without a forward primer... for a reverse primer
cutadapt -a TYRATCAAGAACGAAAGT -O 3 --discard-untrimmed -m 10 -o Test01.assembled.clipped.regR.fasta discard_regF.fasta>> Test01.filter.log

#combine
cat Test01.assembled.clipped.regF.fasta Test01.assembled.clipped.regR.fasta >> Test01.assembled.clipped.fasta
```

## Step 5 - Clean up directory
Move excess files to 'split' directories which were created during the split library step
```
mv Test01*reg* split_Test01/
mv Test01.filter.log split_Test01/
```

## Step 6 - Length filter
Here remove sequences shorter than 150 bp and longer than 500 bps, and rename.
Usage: 
seqlength_cutoff.pl input.fasta min max output.fasta

```
./seqlength_cutoff.pl Test01.assembled.clipped.fasta 150 500 Test01.assembled.clipped.len.fasta
```


## Step 7 - Chimera check 
Remove chimeras using vsearch (uchime) and the PR2 database
Get vsearch here: https://github.com/torognes/vsearch
You will need to acquire reference database that is best for your sample type. I am using the PR2 database. Alternatively, you can change the command to run vsearch chimera checking de novo. This will take longer.

```
vsearch --uchime_ref Test01.assembled.clipped.len.fasta --db /galadriel/sarah/PR2/pr2.qiime.fasta --uchimeout Test01.uchimeinfo_ref --chimeras Test01.chimeras_ref.fasta --strand plus --nonchimeras Test01.assembled.clipped.len.nc.final.fasta 

#Move excess files (chimeras) to split directories
mv Test01.chimeras_ref.fasta split_Test01/
mv Test01.uchimeinfo_ref split_Test01/

```

## Step 8

Combine all reads together. Next run through will add the next set of sequences.
 
```
cat Test01.assembled.clipped.len.nc.final.fasta >> allseqs_test.fasta
```

## Step 9 (optional)
If you have installed abyss or have another way to get quick stats on fasta files, use it here!
It is important to track how many sequences were lost during the QC process. Run this bit of code to obtain separate .txt files for each step.
```
abyss-fac -v *_merged.fastq >>stats_merged.txt
abyss-fac -v *.merged.Q30.fasta >> stats_Q30.txt
abyss-fac -v *.assembled.clipped.fasta >> stats_primerclipped.txt
abyss-fac -v *.assembled.clipped.len.fasta >> stats_length_cutoff.txt
abyss-fac -v *.assembled.clipped.len.nc.final.fasta >> stats_chimeras.txt
```
Text file outputs from this can then be compiled using an R script to generate a QC stats file.
```
#Start R
stats<-c("stats_merged.txt", "stats_Q30.txt", "stats_primerclipped.txt", "stats_length_cutoff.txt", "stats_chimeras.txt")

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

### Next steps - OTU clustering
I've included the test file result "allseqs_test.fasta" here as well.

# Considerations for OTU clustering
Depending on the questions you want to ask your data, decide about what type of OTU clustering you want to do (both algorithm and percent similarity).
 
QIIME is a great resource for tutorials on OTU clustering - http://qiime.org/tutorials/otu_picking.html#running-the-otu-picking-workflows

In this example I will use open reference OTU picking (Rideout et al. 2014, PeerJ), via QIIME and the PR2 database.
```
#OTU clustering with uclust, will make a new directory (pick_open)
pick_open_reference_otus.py -i allseqs_test.fasta -o pick_open -m uclust -r /galadriel/sarah/PR2/pr2.qiime.fasta --suppress_step4 --suppress_taxonomy_assignment
cd pick_open
#assign taxonomy using PR2 database, creates a new directory (uclust_taxonomy)
assign_taxonomy.py -i rep_set.fna -t /galadriel/sarah/PR2/ids.names.2.txt -r /galadriel/sarah/PR2/pr2.qiime.fasta -m uclust -o uclust_taxonomy
#make an OTU table, and convert to txt
make_otu_table.py -i final_otu_map.txt -o V4_tagseq_test.biom -t uclust_taxonomy/rep_set_tax_assignments.txt 
biom convert -i V4_tagseq_test.biom -o V4_tagseq_test.txt --to-tsv --header-key=taxonomy
```
### Use a text editor to remove the comment from the second line
Keep the "Constructed from biom file" commented out line.
See other repository: "PreliminaryFigures_V4_tagseq" for generating preliminary figures from OTU tables.

## Contributing
With helpful guidance and input from Jay Liu

## Citations
QIIME, Caporaso et al. 2010. Nature Methods; doi:10.1038/nmeth.f.303
cutadapt, Martin, 2017; DOI:10.14806/ej.17.1.200 

## Version/update notes
Initial upload 06/2017



