#As of June 2017
This project is under construction! It is a work in progress. If you would like more information, contact me directly at sarah.hu@usc.edu
Also, please feel free to contact me with suggestions or feedback on this pending repo
#V4 tag sequence QC automation
Protocol describes quality checking process for raw fastq sequences. Includes some discussion for alternate approaches. End product is ready for OTU clustering.

I've also included suggestions for OTU clustering at the end (specific for 18S data).

##Getting started
Step by step instructions below. Included perl script here automates the entire process. Directions on automating are included here (protocols.io) link. For first timers, we suggest going through all steps below with one set of fastq samples (e.g. R1 and R2s) to familiarize yourself with process. 

This protocol (and my preferred method) is to perform as much quality checking on each individual sample before combining for downstream OTU clustering (or other analysis). This way, I can track any samples that may be not ideal for downstream analysis.

To follow along with below step by step instructions:
- Download and unzip test_fastq_files.zip
- Download other scripts: create.map.pl, seqlength_cutoff.pl
- Download PR2 database from here < need to figure out how to get pr2 db

###Prerequisites
Programs required:

QIIME - I'm using v.1.9.1
Install here: http://qiime.org/install/install.html
Make sure fastqjoin is installed. Alternatively you can use PEAR to merger PE sequences: https://sco.h-its.org/exelixis/web/software/pear/doc.html#installing

cutadapt - I'm using v.1.10: http://cutadapt.readthedocs.io/en/stable/guide.html

vsearch - v1.11.1: https://github.com/torognes/vsearch

To know about your specific samples:
- Forward and reverse primer sequences, I'm using V4 Stoeck et al. 2010 primers 
- Preferred reference database, here, I'm using the PR2 reference database

##Step 1 - Create map file
Map files in QIIME serve to de-multiplex reads (based on barcodes and indices) and are used in 'split_library' QIIME step. Here, I've split up the V4 forward primer to fulfill the necessary barcode / index space required in the QIIME mapping file. However, we aren't using this feature of the mapping file. 
 
Use this QIIME command to check mapping file:
http://qiime.org/1.6.0/scripts/check_id_map.html

```
perl create.pl
```
This reads the prefix.txt file to create map files specific to the test fastq files for this example.

Following steps only refer to "Test01" sample.


##Step 2 - Merge paired end (R1 and R2) reads with a 20 bp overlap

```
join_paired_ends.py -f Test01_L001_R1_001.fastq -r Test01_L001_R2_001.fastq -o joined_Test01 -j 20

mv joined_Test01/fastqjoin.join.fastq Test01_merged.fastq
```

##Step 3 - Initial Q score filter
Perform split library command in QIIME, filter sequences with Q score 30

```
split_libraries_fastq.py -i Test01_merged.fastq -m Test01_map.txt --barcode_type 'not-barcoded' --sample_ids Test01 -q 29 -n 0 -o split_Test01
mv split_Test01/seqs.fna Test01.merged.Q30.fasta
```

##Step 4 - Remove primers

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

###Step 5 - Clean up directory
Move excess files to 'split' directories which were created during the split library step
```
mv Test01*reg* split_Test01/
mv Test01.filter.log split_Test01/
```

###Step 6 - Length filter
Here remove sequences shorter than 150 bp and longer than 500 bps, and rename.
Usage: 
seqlength_cutoff.pl input.fasta min max output.fasta

```
seqlength_cutoff.pl Test01.assembled.clipped.fasta 150 500 Test01.assembled.clipped.len.fasta
```


###Step 7 - Chimera check 
Remove chimeras using vsearch (uchime) and the PR2 database
Get vsearch here: https://github.com/torognes/vsearch
You will need to acquire reference database that is best for your sample type. I am using the PR2 database. Alternatively, you can change the command to run vsearch chimera checking de novo. This will take longer.

```
vsearch --uchime_ref Test01.assembled.clipped.len.fasta --db /galadriel/sarah/PR2/pr2.qiime.fasta --uchimeout Test01.uchimeinfo_ref --chimeras Test01.chimeras_ref.fasta --strand plus --nonchimeras Test01.assembled.clipped.len.nc.final.fasta 

#Move excess files (chimeras) to split directories
mv Test01.chimeras_ref.fasta split_Test01/
mv Test01.uchimeinfo_ref split_Test01/

```

###Step 8

Combine all reads together. Next run through will add the next set of sequences.
 
```
cat Test01.assembled.clipped.len.nc.final.fasta >> allseqs_test.fasta
```

## Contributing
With helpful guidance and input from Jay Liu

## Citations
QIIME, Caporaso et al. 2010. Nature Methods; doi:10.1038/nmeth.f.303
cutadapt, Martin, 2017; DOI:10.14806/ej.17.1.200 

## Version/update notes
Initial upload 06/2017



