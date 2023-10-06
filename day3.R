library(tidyverse) 
library(reshape2)
library(vegan)
library(rioja)
library(ggplot2)
library(dplyr)
library(gghighlight)

#First, set working directory

setwd("~/course/wdir/mapping/plots/")

# Import data

df <- read_csv("~/course/wdir/mapping/metaDMGresults.csv")

# Count number of columns in dataframe, just to familiarize your self with it.

ncol(df)

# Print all sample names

unique(df$sample)

# Import metadata

metaDATA <- read.delim("~/course/data/shared/metadata/metadata.tsv")

# Replace header "sample_name" with "sample"

colnames(metaDATA)[colnames(metaDATA) == "sample_name"] <- "sample"
colnames(metaDATA)[colnames(metaDATA) == "years_bp"] <- "YearsBP"

# Merge the metadata and dataframe by sample

dt <- merge(df, metaDATA, by = "sample")

# Now check that new columns have been added to the dt dataframe

ncol(df) < ncol(dt)

# Print all different ages

unique(dt$YearsBP)

# Print all column names of the table

colnames(dt)

# Maximum amount of damage (filtered for)

DamMin2 = 0.00

# MAP_Significance

MapSig2 = 0

# Minimum reads for parsing taxa

MinRead2 = 100

# Minimum mean read length

MinLength = 35

# Subsetting the table by plant (Viridiplantae) genus using grepl and filter, parameters you need to set and possible add more?? Please think about how your group think it should be done. But first would probably be best to explore the data with less stringent filters, and plot these, then later on the basis of the data make decisions on how you want it filtered.

dt2 <- dt %>% filter(MAP_damage > DamMin2, N_reads >= MinRead2, mean_L > MinLength, MAP_significance  > MapSig2,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", sample))

# Now plot your data, you can add other variables in the gghighlight to illustrate different cut-offs, or other types of values etc.

pdf(file = "aeCourse.DNAdamageModelJitterPlot.pdf", width = 8, height = 4)
ggplot() +
  geom_jitter(data = dt2, aes(x=as.numeric(YearsBP), y=MAP_damage, size = N_reads), alpha =0.5) +
  gghighlight(N_reads > 500) +
  xlab("Years BP") +
  ylab("DNA damage") +
  labs(color = "Values for taxa with \n>500 reads", size = "Number of reads")
dev.off()

# Plot plant taxa, highlight taxa with more than 500 reads and add the min, max and median. save as pdf

pdf(file = "aeCourse.DNAdamageLRJitterPlot.pdf", width = 8, height = 4)
ggplot() +
  geom_jitter(data = dt2, aes(x=as.numeric(YearsBP), y=MAP_damage, size = MAP_significance), alpha =0.5) +
  gghighlight(N_reads > 500) +
  xlab("Years BP")+
  ylab("DNA damage") +
  labs(color = "Values for taxa with \n>500 reads", size = "Significance \nfor Taxa with >500 reads")
dev.off()

# Create filtered table for DNA damage model

filtered_data <- dt2 %>% filter(N_reads >= 500)

# Now you should make decisions on how the data should be filtered. You can use the

MapSig3 = 3
MinRead3 = 100
MinLength3 = 35

# Subsetting the table using grepl and filter, parameters you need to set and possible add more?

filtered_data_metazoan <- dt %>% filter(N_reads >= 10, mean_L > MinLength3, MAP_significance  > MapSig3,  grepl("Metazoa",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", YearsBP))
unfiltered_data_metazoan <- dt %>% filter(N_reads >= 2, mean_L > MinLength3, MAP_significance  > 1,  grepl("Metazoa",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", YearsBP))
filtered_data_viridiplantae <- filtered_data %>% filter(N_reads >= 100, mean_L > MinLength3, MAP_significance  > MapSig3,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", YearsBP))
unique(filtered_data_viridiplantae$YearsBP)

# You can also make a smaller table with only values of your choice example below.

select(filtered_data_viridiplantae, tax_name, MAP_damage, MAP_significance, N_reads, YearsBP)

# Count number of unique plant taxa

unique(filtered_data_viridiplantae$tax_name)

# Make wide table for downstream plot and data wrangling of the plants.

data_wide_plants <- dcast(filtered_data_viridiplantae, tax_name ~ YearsBP, value.var="N_reads", fun.aggregate = sum)
n <- ncol(data_wide_plants)
b2 <- data_wide_plants[,2:n]
rownames(b2) <- data_wide_plants$tax_name
b2[is.na(b2)] <- 0 #set all NAs as zeros

# Prints sum of samples and taxa and depths names

colSums(b2)
rowSums(b2)
colnames(b2)

# Test, if this one fails it might have text in the number of reads coloumn

b2[is.na(b2)]=0

# Create percentage table, by taking the number of columns and divides these with sum of all reads in the column, remove the -1 and you will plot the Nreads column on later plots

i=ncol(b2)
b3=as.matrix(b2[,seq(1,i)])  

b4 <- prop.table(data.matrix(b3), margin=2)*100 # makes proportion table, needs 2 margins e.g. header and 1st row names
colSums(prop.table(b4, margin=2)*100) # prints sum of column, should give 100 for each

# Next we will transpose the table, for plotting it as a strat.plot

b5 <- t(b4)

# and set the variable z to be the years BP which is now row headers.

z <- as.numeric(rownames(b5)) # depth/depth

# and plot it on a stratigraphic plot (typical pollen type plot)

pdf(file = "aeCourse.Stratplot_Plants_area.pdf", width = 15, height = 5)
strat.plot(b5, y.rev=TRUE, plot.line=TRUE, plot.poly=TRUE, plot.bar=FALSE, lwd.bar=10, sep.bar=TRUE, scale.percent=TRUE, xSpace=0.01, x.pc.lab=TRUE, x.pc.omit0=TRUE, srt.xlabel=45, las=2, exag=TRUE, exag.mult=5, ylabel = "years BP",  yvar = z)

dev.off()

# Now we will convert the wide table into a long table format to plot it with ggplot

y <- ncol(b5)
b6 <- melt(b5[,1:y])
sapply(b6, class)
colnames(b6) <- c("YearsBP","Taxa", "percentage")

p1 <- ggplot(b6, aes(y=Taxa, x=YearsBP, fill=percentage)) +   geom_tile(colour="lightgrey") +
  theme_minimal() + scale_fill_gradient(low="white", high="darkgreen") + scale_y_discrete(limits=rev)
p1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust =1)) + ggtitle("percentage of taxa plotted as heatmap") +
  xlab("YearsBP") + ylab("YearsBP") + labs(fill = "percentage %")
