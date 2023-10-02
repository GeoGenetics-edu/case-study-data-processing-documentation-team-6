### R script
library(ggplot)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
read_len <- read_tsv(args[1], header = F, colnames = c("count", "length"))
read_len %>% 
ggplot(.,aes(x = count)) + 
    geom_histogram()

ggsave(args[2])
