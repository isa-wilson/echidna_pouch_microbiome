# Identification of contaminants in echidna pouch microbiome dataset 
# Isabella Wilson 2024-12-20
# Script adapted from Raphael Eisenhofer 2020-08-04

library(tidyverse)
library(readr)
library(phyloseq)
library(ggplot2)
library(decontam)
library(scales)
library(qiime2R)
library(svglite)
library(cowplot)

# Note: before starting, edit .tsv metadata file: 
#       - change "sample-id" column name to "SampleID" 
#       - add "isControl" column with "TRUE" or "FALSE" as values

# Import metadata
metadata<-read_tsv("sample-metadata.tsv")

Table <- read_qza("freq-filtered-table.qza") # Note - this table has had singletons removed

# Import silva database and phylogenetic tree 

silva_taxonomy <- read_qza("taxonomy.qza")

rooted_tree <- read_qza("rooted-tree.qza")

physeq <- phyloseq(otu_table(Table$data, taxa_are_rows = TRUE))

taxtable<-silva_taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) 

SBEPhy<-phyloseq(otu_table(Table$data, taxa_are_rows = T), phy_tree(rooted_tree$data), tax_table(as.data.frame(taxtable) %>% column_to_rownames("Feature.ID") %>% as.matrix()), sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sample-id")))

# Generate library size plot:
SBE_df <- as.data.frame(sample_data(SBEPhy))

SBE_df$LibrarySize<-sample_sums(SBEPhy)

SBE_df <- SBE_df[order(SBE_df$LibrarySize),]
SBE_df$Index <- seq(nrow(SBE_df))

ggplot(data=SBE_df, aes(x=Index, y=LibrarySize, color="SampleType")) +geom_point(size=1)
ggsave(filename = "decontam/LibrarySize.svg", dpi = 300)


# Perform prevalence-analysis
sample_data(SBEPhy)$is.neg <- sample_data(SBEPhy)$isControl == "TRUE"

SBE_contamdf.prev <- isContaminant(SBEPhy, method="prevalence", neg="is.neg")

SBE_contamdf.prev05 <- isContaminant(SBEPhy, method="prevalence", neg="is.neg", threshold=0.5)

table(SBE_contamdf.prev$contaminant)
table(SBE_contamdf.prev05$contaminant)

# Generate prev-prev plot
SBE_Phy.pa <- transform_sample_counts(SBEPhy, function(abund) 1*(abund>0))

SBE_Phy.pa.controls <- prune_samples(sample_data(SBE_Phy.pa)$isControl == "TRUE", SBE_Phy.pa)

SBE_Phy.pa.samples <- prune_samples(sample_data(SBE_Phy.pa)$isControl == "FALSE", SBE_Phy.pa)

SBEPhy_df.pa <- data.frame(pa.pos=taxa_sums(SBE_Phy.pa.samples), pa.neg=taxa_sums(SBE_Phy.pa.controls), contaminant=SBE_contamdf.prev$contaminant)

ggplot(data=SBEPhy_df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point(size=1) + geom_jitter() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (Samples)") 
ggsave(filename = "decontam/decontam-prev-prev.svg", dpi = 300)

# Generate decontam histogram
SBEPhy_df05.pa <- data.frame(pa.pos=taxa_sums(SBE_Phy.pa.samples), pa.neg=taxa_sums(SBE_Phy.pa.controls), contaminant=SBE_contamdf.prev05$contaminant)

ggplot(data=SBEPhy_df05.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + xlab("Prevalence (Negative Controls)") + ylab("Prevalence (Samples)") 

ggplot(data = SBE_contamdf.prev, aes(x=p)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = 'decontam Score', y='Number of species')

ggsave(filename = "decontam/decontam-score-histo.svg", dpi = 300)

# Get taxonomy of putative contaminants
tax <- as(tax_table(SBEPhy), "matrix")

# Coerce into a dataframe
tax.df <- data.frame(tax)

# Merge dataframes to obtain taxonomic information
test.merge <- merge(SBE_contamdf.prev05, tax.df, by=0, sort=FALSE)

# Generate two subsets of taxa: contaminants and non-contaminants
contaminant.ASVs <- subset(test.merge, contaminant == "TRUE")

non.contaminant.ASVs <- subset(test.merge, contaminant == "FALSE")

# Export as TSVs for table filtering in Qiime2 
write.table(test.merge, "decontam/decontam.results.0.5.ASVs.txt", sep="\t")
write.table(contaminant.ASVs, "decontam/decontam.0.5.contaminant.ASVs.txt", sep="\t")
write.table(non.contaminant.ASVs, "decontam/decontam.0.5.non.contaminant.ASVs.txt", sep="\t")
