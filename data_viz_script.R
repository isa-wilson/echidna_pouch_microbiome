### Scripts for microbiome plots in R ###

# Isabella Wilson 17/11/23

library(tidyverse)
library(qiime2R)
library(phyloseq)
library(fantaxtic)
library(microViz)
library(microbiome)

# Read in Qiime files
metadata <- read_tsv("sample-metadata-pouch.tsv")

# Convert to dataframe
metadata.df <- as.data.frame(metadata)

table <- read_qza("filtered-table-pouch.qza")$data

taxonomy<-read_qza("taxonomy.qza")$data %>% parse_taxonomy()

# Convert to matrix
taxonomy.matrix <- as.matrix(taxonomy)

# Create phyloseq object
OTU = otu_table(table, taxa_are_rows=TRUE)
TAX = tax_table(taxonomy.matrix)
sampledata = sample_data(metadata.df)
rownames(sampledata) <- metadata.df$`sample-id`

physeq = phyloseq(OTU,TAX,sampledata)

# Rename "NA"/"uncultured" taxa with last identified rank
physeq.fix <- tax_fix(physeq.fix, unknowns = c("Ambiguous_taxa", "uncultured", "NA"))
 
# Top 15 genera bar plot                     
physeq.fix %>%
  comp_barplot("Genus", n_taxa = 15, merge_other = FALSE, bar_outline_colour = "#f9f9f9") +
  facet_grid(~BreedingStatus, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45) + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black',))
  )

### Fancy nested plots ###

# Get the most abundant phyla and the most abundant families within those phyla
top_nested <- nested_top_taxa(physeq,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Genus",
                              n_top_taxa = 4, 
                              n_nested_taxa = 3)

# Plot the relative abundances at two levels
basic.plot <- plot_nested_bar(top_nested$ps_obj,
                              top_level = "Phylum",
                              nested_level = "Genus",
                              nested_merged_label = "NA and other <tax>",
                              legend_title = "Phylum and Genus") + guides(fill = guide_legend(ncol = 1))


# Add faceting
basic.plot + facet_grid(~BreedingStatus,
                        scales = "free", space = "free") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 8, 
                                  face = "bold"),
        legend.key.size = unit(10, 
                               "points")) +
  theme(axis.text.x = element_text(angle = 90)) + 
  theme(strip.background =element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black',))



### PCOA plots ###

# Read in data

table.nl <- read_qza("filtered-table-pouch-nonlactating.qza")$data
taxonomy.nl <- read_qza("taxonomy.qza")$data %>% parse_taxonomy()

# Convert to matrix
taxonomy.nl.matrix <- as.matrix(taxonomy.nl)
rooted.tree <- read_tree("export/tree.nwk")
metadata.nl <- read_tsv("sample-metadata-pouch-nonlactating.tsv")

# Convert to dataframe
metadata.nl.df <- as.data.frame(metadata.nl)

# Create physeq object
OTU.nl = otu_table(table.nl, taxa_are_rows=TRUE)
TAX.nl = tax_table(taxonomy.nl.matrix)
sampledata.nl = sample_data(metadata.nl.df)
rownames(sampledata.nl) <- metadata.nl.df$`sample-id`

physeq.nl = phyloseq(OTU.nl,TAX.nl,sampledata.nl,rooted.tree)

# Plot PCoA
ordu = ordinate(physeq.nl, "PCoA", "unifrac", weighted=TRUE)
PCOA.1.2 = plot_ordination(physeq.nl, ordu, color="CaptivityStatus", axes = c(1,2))
PCOA.1.3 = plot_ordination(physeq.nl, ordu, color="CaptivityStatus", axes = c(1,3))

# Axes 1 and 2
PCOA.1.2 +
  scale_colour_manual(values = c("Captive" = "#f27304", "Wild" = "#008000")) +
  geom_point(size = 3.5) +
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black", linewidth = 0.75),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
  ) 
ggsave(filename = "Rplot-captivity-PCoA-axes-1-2.svg", dpi = 300)

# Axes 1 and 3
PCOA.1.3 +
scale_colour_manual(values = c("Captive" = "#f27304", "Wild" = "#008000")) +
geom_point(size = 3.5) +
theme_bw() +
theme(
  panel.border = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  axis.line = element_line(colour = "black", linewidth = 0.75),
  axis.ticks = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  legend.position = "none",
) 
ggsave(filename = "Rplot-captivity-PCoA-axes-1-3.svg", dpi = 300)
