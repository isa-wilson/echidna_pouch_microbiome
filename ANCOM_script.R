### Differential abundance with ANCOM-BC ###

library(ANCOMBC)
library(tidyverse)
library(qiime2R)
library(phyloseq)

# import data

metadata <- read_tsv("sample-metadata-pouch.tsv")
metadata.df <- as.data.frame(metadata)

table <- read_qza("filtered-table-pouch.qza")$data

taxonomy<-read_qza("taxonomy.qza")$data %>% parse_taxonomy()
taxonomy.matrix <- as.matrix(taxonomy)

# create phyloseq object
OTU = otu_table(table, taxa_are_rows=TRUE)
TAX = tax_table(taxonomy.matrix)
sampledata = sample_data(metadata.df)
rownames(sampledata) <- metadata.df$`sample-id`

physeq = phyloseq(OTU,TAX,sampledata)

# reduce phyloseq object to taxonomic rank of interest
physeq_genus <- tax_glom(physeq, taxrank = "Genus")

physeq_phylum <- tax_glom(physeq, taxrank = "Phylum")

# run ANCOM
out.genus = ancombc(
  phyloseq = physeq_genus, 
  formula = "SampleType", 
  p_adj_method = "holm", 
  zero_cut = 0.90, # by default prevalence filter of 10% is applied
  lib_cut = 0, 
  group = "SampleType", 
  struc_zero = FALSE, 
  neg_lb = FALSE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)

out.phylum = ancombc(
  phyloseq = physeq_phylum, 
  formula = "SampleType", 
  p_adj_method = "holm", 
  zero_cut = 0.90, # by default prevalence filter of 10% is applied
  lib_cut = 0, 
  group = "SampleType", 
  struc_zero = FALSE, 
  neg_lb = FALSE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)

res = out.genus$res

