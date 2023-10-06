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
- step 7: Visualize the read length distribution

### The shell script
```
#!/bin/bash -l

# step 1: load the working environment
conda activate day1

# MAIN LOOP: Executing for each file
# TEST ONE FILE: for f in ~/course/data/day2/fastq/PRI-TJPGK-CATN-96*
for f in ~/course/data/day2/fastq/PRI-TJPGK-*
do

# step 2: Make symbolic link
ln -s $f .
BASENAME=$(basename $f)

# step 3: Remove duplicates and reads less than 30 bp
vsearch --fastx_uniques $BASENAME --fastqout ./${BASENAME/.fq.gz/.vs.fq} --minseqlength 30 --strand both

# step 4: Compress output
gzip ${BASENAME/.fq.gz/.vs.fq}

# step 5a: Align reads to aegenmoics.db and covert to bam format
bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U ${BASENAME/.fq.gz/.vs.fq.gz} --no-unal | samtools view -bS - > ${BASENAME/.fq.gz/.bam}

# step 5b: Sort the bam file by read name
samtools sort -n ${BASENAME/.fq.gz/.bam} > ${BASENAME/.fq.gz/.sort.bam}

# step 6: Calculate read lengths
zcat ${BASENAME/.fq.gz/.vs.fq.gz} | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l"\t"lengths[l]}}' | sort -n > ${BASENAME/.fq.gz/.read_lengths.txt}

# step 7: visualization 
Rscript histogram.R ${BASENAME/.fq.gz/.read_lengths.txt} ${BASENAME/.fq.gz/.read_lengths.png}

done
```
### The Rscript (step7)
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


## Assignment 2

Compute DNA damage and statistics for all 5 samples
```
conda activate metaDMG

metaDMG config *.sort.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp

# Edit config to change parameters for exploration
vim config.yaml

metaDMG compute config.yaml

# And output as a csv
metaDMG convert --output metaDMGresults.csv --add-fit-predictions
```

### Testing similarity score parameters
We ran metaDMG with minimum similarity scores of 0.8 and 0.95. Somewhat counterintuitively, less taxa were present in the resulting data with a similarity score 0.8 compared to 0.95. We hypothesize that the lower similarity score minimum allows for more alignments which can lead to the LCA assigning the reads to a higher taxonomic level. 

### Testing bayesian estimation
... still running ...

## Assignment 3

### aDNA authentication 

Answers here..

### euka
Alternative assignment from holi pipeline. 

`
~/course/data/vgan/bin/vgan euka -fq1 <(zcat *vs.fq.gz) -o all_samples -t 5 --euka_dir euka_dir/
`
Ursidae was not detected because the reads added to the query dataset were extracted from a competitively matched dataset so conserved regions of the mitochondria were missing, causing `euka` to label the poor breadth of coveraged reads as contaminants.

Found Bovidae which was not present in holi results. This is due to the absence of Bovidae in the reference database used in holi:

`zgrep NC_001941 ~/course/data/shared/taxonomy/acc2taxid.map.gz`

### pathPhynder
```
conda activate day2

mamba install unzip
unzip pathphynder_analysis.zip

snp-sites -v -c -o Aln_mafft_taxa_references.vcf All_ref_only_realigned.fa

mamba install r-stringr

Rscript ./fix_vcf.R Aln_mafft_taxa_references.vcf Aln_mafft_taxa_references_output.vcf

awk '{ if ($1 == "1") $1="consensus";}1' Aln_mafft_taxa_references_output.vcf | sed 's/ /\t/g' > Aln_mafft_taxa_references_fixed_output.vcf

mamba install biopython

python ./get_consensus.txt All_ref_only_realigned.fa Cons_Aln_mafft_taxa_references.fa

conda activate day1
bwa index Cons_Aln_mafft_taxa_references.fa

bwa aln -l 1024 -n 0.001 -t 5 ./Cons_Aln_mafft_taxa_references.fa ./PRI-DJYLO-VM-17-3-1-8-iPCR2_S4_L001_R1_001.lca.txt.Ovis.fq | bwa samse ./Cons_Aln_mafft_taxa_references.fa  - ./PRI-DJYLO-VM-17-3-1-8-iPCR2_S4_L001_R1_001.lca.txt.Ovis.fq | samtools view -F 4 -q 25 -@ 5 -uS - | samtools sort -@ 5 -o PRI-DJYLO-VM-17-3-1-8-iPCR2_S4_L001_R1_001.lca.txt.Ovis.sort.bam

samtools view -c PRI-DJYLO-VM-17-3-1-8-iPCR2_S4_L001_R1_001.lca.txt.Ovis.sort.bam

git clone https://github.com/ruidlpm/pathPhynder.git
touch ~/.bash_profile
vi ~/.bash_profile

source ~/.bash_profile

conda activate r
mamba install r-phytools
mamba install r-optparse

pathPhynder -B -o ./branches.snp ./Mafft_All_BEAST4.nwk ./Mafft_All_references_fixed_consensus.vcf

mamba install -c bioconda samtools
pathPhynder -s prepare -i Mafft_All_BEAST4.nwk -p taxa_pathphynder_tree -f branches.snp -r Cons_Aln_mafft_taxa_references.fa
pathPhynder -s all -t 100 -m transversions -i Mafft_All_BEAST4.nwk -p tree_data/taxa_pathphynder_tree -r Cons_Aln_mafft_taxa_references.fa -b PRI-DJYLO-VM-17-3-1-8-iPCR2_S4_L001_R1_001.lca.txt.Ovis.sort.bam
pathPhynder -s all -t 100 -i Mafft_All_BEAST4.nwk -p tree_data/taxa_pathphynder_tree -r Cons_Aln_mafft_taxa_references.fa -b PRI-DJYLO-VM-17-3-1-8-iPCR2_S4_L001_R1_001.lca.txt.Ovis.sort.bam 

```
### Populations genomics

#### Including polar bear in calculting PC space
- Which populations are separated on the first few PCs?
  
  Polar bears are separated from all other bears on PC1 which dominates the space.
- Where do the ancient bear samples fall, and how can we interpret their position?
  
  The ancient bears samples cluster with all black bears on PC1, Kenai samples on PC2, and East and West on PC4
- Is there any difference in PCA positions between the three bear samples, and if so, what could be the interpretation?
  
  They fall roughly in the same space with UE1605 slighly closer to polar bears on PC1

#### Excluding polar bear in calculating PC space
- Which populations are separated on the first few PCs?
  
  All populations are well differentiated with PC1 and PC2.
- Where do the ancient bear samples fall, and how can we interpret their position?
- Is there any difference in PCA positions between the three bear samples, and if so, what could be the interpretation?

##### Neighborhood joining tree

![fst_tree png-1](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team-6/assets/2534009/afbd67d0-8f16-4303-b85a-2a90261c34dd)


## Day 5
### Pathogen
#### Read mapping
- What is the average coverage we obtained on the main chromosome?​​
- What fraction of the mapped reads were duplicates?
- How much was the coverage reduced?
- Does the coverage look even across all contigs?​​
#### Damage estimates
- Are the damage patterns consistent with ancient Yersinia DNA?​​
  
#### Detailed coverage statistics

