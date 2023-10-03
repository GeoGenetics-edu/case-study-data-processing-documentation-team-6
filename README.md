[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)
# Assignment 
## Assignment 1
### The outline of the script - assignment1.sh
- step 1: Load the working environment
- step 2: Make symbolic link for all fastq file
- step 3: Remove duplicates and reads' length less than 30 bp
- step 4: Compress the output from step 4
- step 5: Align reads to aegenmoics.db, covert alignment output into bam format and sort the bam file
- step 6: Calculate read lengths
- setp 7: Visualize the read length distribution

### The shell script
```
#!/bin/bash -l

# load the working environment
conda activate day1

# MAIN LOOP: Executing for each file
# TEST ONE FILE: for f in ~/course/data/day2/fastq/PRI-TJPGK-CATN-96*
for f in ~/course/data/day2/fastq/PRI-TJPGK-*
do

# Make symbolic link
ln -s $f .
BASENAME=$(basename $f)

# Remove duplicates and reads less than 30 bp
vsearch --fastx_uniques $BASENAME --fastqout ./${BASENAME/.fq.gz/.vs.fq} --minseqlength 30 --strand both

# Compress output
gzip ${BASENAME/.fq.gz/.vs.fq}

# Align reads to aegenmoics.db and covert to bam format
bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U ${BASENAME/.fq.gz/.vs.fq.gz} --no-unal | samtools view -bS - > ${BASENAME/.fq.gz/.bam}

# Sort the bam file by read name
samtools sort -n ${BASENAME/.fq.gz/.bam} > ${BASENAME/.fq.gz/.sort.bam}

# Calculate read lengths
zcat ${BASENAME/.fq.gz/.vs.fq.gz} | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l"\t"lengths[l]}}' | sort -n > ${BASENAME/.fq.gz/.read_lengths.txt}

# visualization 
Rscript histogram.R ${BASENAME/.fq.gz/.read_lengths.txt} ${BASENAME/.fq.gz/.read_lengths.png}

done
```
### The Rscript
```
### R script

#### load R packages
library(ggplot2)

### read the parameters
args = commandArgs(trailingOnly=TRUE)

### load the input file
read_len <- read.table(file = args[1], sep = '\t', header = FALSE, col.names = c("length", "count"))

### visualization
ggplot(read_len, aes(x = length, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Read Length", y = "Count") +
  theme_bw()

### save the figure
ggsave(args[2])
```

### Run the shell script
```
bash assignment1.sh
```
