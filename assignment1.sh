#!/bin/bash -l

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

# Calculate read lengths
zcat ${BASENAME/.fq.gz/.vs.fq.gz} | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l"\t"lengths[l]}}' | sort -n > ${BASENAME/.fq.gz/.read_lengths.txt}

# visualization 
Rscript histogram.R ${BASENAME/.fq.gz/.read_lengths.txt} ${BASENAME/.fq.gz/.read_lengths.png}

done
