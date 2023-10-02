### R script

#### load R packages
library(ggplot2)
library(tidyverse)

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
