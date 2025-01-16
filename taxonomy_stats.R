# Statistical analysis of echidna pouch microbiome taxonomy data
# Isabella Wilson 2024-12-20

library(tidyverse)
library(readr)
library(qiime2R)
library(phyloseq)
library(microViz)

# Read in data ----  
metadata <- read_tsv("sample-metadata-pouch.tsv")
metadata.df <- as.data.frame(metadata) # Convert to dataframe
table <- read_qza("filtered-table-pouch.qza")$data
taxonomy<-read_qza("taxonomy.qza")$data %>% parse_taxonomy()
taxonomy.matrix <- as.matrix(taxonomy) # Convert to matrix

# Create separate table and metadata objects according to lactation status
table.lact <- table[, c("MSSBE114","MSSBE115","MSSBE116","MSSBE117","MSSBE118","MSSBE119")]
table.nonlact <- table[, c("3B","2P","NBSBE3","6B","2B","9B","10B","7B","1B","NBSBE2","7PM","5B","IWSBE5","NBSBE1","8B","MSSBE113")]

metadata.sort <- arrange(metadata, desc(metadata$SampleType))
metadata.lact <- metadata.sort[17:22,]
metadata.nonlact <- metadata.sort[1:16,]

# Create phyloseq objects ----
OTU.lact = otu_table(table.lact, taxa_are_rows=TRUE)
TAX = tax_table(taxonomy.matrix)
sampledata.lact = sample_data(metadata.lact)
rownames(sampledata.lact) <- metadata.lact$`sample-id`

physeq.lact = phyloseq(OTU.lact,TAX,sampledata.lact)

OTU.nonlact = otu_table(table.nonlact, taxa_are_rows=TRUE)
TAX = tax_table(taxonomy.matrix)
sampledata.nonlact = sample_data(metadata.nonlact)
rownames(sampledata.nonlact) <- metadata.nonlact$`sample-id`

physeq.nonlact = phyloseq(OTU.nonlact,TAX,sampledata.nonlact)

# Calculate averages at desired tax rank ----
## Phylum ----
### lactating ----
phylum.lact = tax_glom(physeq.lact, "Phylum", NArm = FALSE) %>% transform_sample_counts(function(x) {x * 100/sum(x)})

# calculate per-row average in the OTU table
df.phy.lact = data.frame(Phylum = tax_table(phylum.lact)[,"Phylum"], Mean = rowMeans(otu_table(phylum.lact)), row.names = NULL)
df.phy.lact = df.phy.lact[order(-df.phy.lact$Mean),]

### non-lactating ----
phylum.nonlact = tax_glom(physeq.nonlact, "Phylum", NArm = FALSE) %>% transform_sample_counts(function(x) {x * 100/sum(x)})

# calculate per-row average in the OTU table
df.phy.nonlact = data.frame(Phylum = tax_table(phylum.nonlact)[,"Phylum"], Mean = rowMeans(otu_table(phylum.nonlact)), row.names = NULL)
df.phy.nonlact = df.phy.nonlact[order(-df.phy.nonlact$Mean),]

head(df.phy.lact)
head(df.phy.nonlact)

## Genus ----
### lactating ----
physeq.lact.fix <- tax_fix(physeq.lact, unknowns = c("Ambiguous_taxa", "uncultured", "NA", "Unknown_Family")) # replaces unclassified genera with more informative name
genus.lact = tax_glom(physeq.lact.fix, "Genus", NArm = FALSE) %>% transform_sample_counts(function(x) {x * 100/sum(x)})

# calculate per-row average in the OTU table
df.gen.lact = data.frame(Genus = tax_table(genus.lact)[,"Genus"], Mean = rowMeans(otu_table(genus.lact)), row.names = NULL)
df.gen.lact = df.gen.lact[order(-df.gen.lact$Mean),]

### non-lactating ----
physeq.nonlact.fix <- tax_fix(physeq.nonlact, unknowns = c("Ambiguous_taxa", "uncultured", "NA", "Unknown_Family"))
genus.nonlact = tax_glom(physeq.nonlact.fix, "Genus", NArm = FALSE) %>% transform_sample_counts(function(x) {x * 100/sum(x)})

# calculate per-row average in the OTU table
df.gen.nonlact = data.frame(Genus = tax_table(genus.nonlact)[,"Genus"], Mean = rowMeans(otu_table(genus.nonlact)), row.names = NULL)
df.gen.nonlact = df.gen.nonlact[order(-df.gen.nonlact$Mean),]

head(df.gen.nonlact)
head(df.gen.lact)
