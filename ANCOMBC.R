# Differential abundance analysis with ANCOMBC
# Isabella Wilson 2024-12-20

library(ANCOMBC)
library(tidyverse)
library(qiime2R)
library(phyloseq)
library(microViz)
library(ggrepel)
library(svglite)

# Prepare data ----

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
  
  # reorder "SampleType" variable in physeq object
  sample_data(physeq)$SampleType <- as.factor(sample_data(physeq)$SampleType)
  sample_data(physeq)$SampleType <- relevel(sample_data(physeq)$SampleType, "Non-lactating pouch")
  
  # reduce phyloseq object to taxonomic rank of interest
  physeq_phylum <- tax_glom(physeq, taxrank = "Phylum")
  physeq_family <- tax_glom(physeq, taxrank = "Family")
  physeq_genus <- tax_glom(physeq, taxrank = "Genus")

# Run ANCOM ----

## Phylum ----
out.phylum = ancombc(
  data = physeq_phylum, 
  formula = "SampleType", 
  p_adj_method = "BH", 
  lib_cut = 0, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = FALSE
)

# convert output to dataframe
res.phylum.df = data.frame(taxon = out.phylum$res$lfc$taxon,
                           Phylum = out.phylum$res$lfc$taxon,
                           lfc = out.phylum$res$lfc$`SampleTypeLactating pouch`, 
                           se = out.phylum$res$se$`SampleTypeLactating pouch`,
                           W = out.phylum$res$W$`SampleTypeLactating pouch`, 
                           p_val = out.phylum$res$p_val$`SampleTypeLactating pouch`, 
                           q_value = out.phylum$res$q_val$`SampleTypeLactating pouch`, 
                           Diff_ab =  out.phylum$res$diff_abn$`SampleTypeLactating pouch`)

# change feature IDs to phylum names
res.phylum.df$Phylum = taxonomy$Phylum[match(res.phylum.df$Phylum ,rownames(taxonomy))]

## Family ----
out.family = ancombc(
  data = physeq_family, 
  formula = "SampleType", 
  p_adj_method = "BH", 
  lib_cut = 0, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = FALSE
)

# convert output to dataframe
res.family.df = data.frame(taxon = out.family$res$lfc$taxon,
                           Family = out.family$res$lfc$taxon,
                           lfc = out.family$res$lfc$`SampleTypeLactating pouch`, 
                           se = out.family$res$se$`SampleTypeLactating pouch`,
                           W = out.family$res$W$`SampleTypeLactating pouch`, 
                           p_val = out.family$res$p_val$`SampleTypeLactating pouch`, 
                           q_value = out.family$res$q_val$`SampleTypeLactating pouch`, 
                           Diff_ab =  out.family$res$diff_abn$`SampleTypeLactating pouch`)

# change feature IDs to family names
res.family.df$Family = taxonomy$Family[match(res.family.df$Family ,rownames(taxonomy))]

## Genus ----

# First rename unclassified genera to next highest taxonomic rank
physeq.fix <- tax_fix(physeq, unknowns = c("Ambiguous_taxa", "uncultured", "NA", "Unknown_Family"))

sample_data(physeq.fix)$SampleType <- as.factor(sample_data(physeq.fix)$SampleType)
sample_data(physeq.fix)$SampleType <- relevel(sample_data(physeq.fix)$SampleType, "Non-lactating pouch")

# Glom to genus level
physeq_genus.fix <- tax_glom(physeq.fix, taxrank = "Genus")

# Run ANCOM
out.genus.fix = ancombc(
  data = physeq_genus.fix, 
  formula = "SampleType", 
  p_adj_method = "BH", 
  lib_cut = 0, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = FALSE
)

# Convert ANCOM output into dataframe
res.genus.df.fix = data.frame(taxon = out.genus.fix$res$lfc$taxon,
                              Genus = out.genus.fix$res$lfc$taxon,
                              lfc = out.genus.fix$res$lfc$`SampleTypeLactating pouch`, 
                              se = out.genus.fix$res$se$`SampleTypeLactating pouch`,
                              W = out.genus.fix$res$W$`SampleTypeLactating pouch`, 
                              p_val = out.genus.fix$res$p_val$`SampleTypeLactating pouch`, 
                              q_value = out.genus.fix$res$q_val$`SampleTypeLactating pouch`, 
                              Diff_ab =  out.genus.fix$res$diff_abn$`SampleTypeLactating pouch`)

# Add genus names to ANCOM output dataframe
tax_table.fix = data.frame(Genus = physeq_genus.fix@tax_table[,"Genus"]) #pull genus names and associated tax IDs from physeq object
res.genus.df.fix$Genus = tax_table.fix$Genus[match(res.genus.df.fix$Genus ,rownames(tax_table.fix))] #convert taxIDs to genus names in ANCOM output dataframe

# Add column describing whether result is associated with lactating/non-lactating pouch
res.genus.df.fix$association = ifelse(res.genus.df.fix$q_value < 0.05 & res.genus.df.fix$lfc > 0, "Lactating", "Not significant")
res.genus.df.fix$association = ifelse(res.genus.df.fix$q_value < 0.05 & res.genus.df.fix$lfc < 0, "Non-lactating",res.genus.df.fix$association)

# Lollipop plot ----
## Phylum ----
res.phylum.sig = res.phylum.df[res.phylum.df$q_value < 0.05,] # separate out significant results

lollipop.phylum = ggplot(res.phylum.sig, aes(lfc, Phylum, color = association ))+
                  geom_point() +theme_bw() +
                  geom_segment(aes(x=0, xend=lfc, y=Phylum, yend=Phylum, color = association)) +
                  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +ylab(NULL)# +

lollipop.phylum

## Family ----
res.family.sig = res.family.df[res.family.df$q_value < 0.05,] # separate out significant results

lollipop.family = ggplot(res.family.sig, aes(lfc, Family, color = association ))+
  geom_point() +theme_bw() +
  geom_segment(aes(x=0, xend=lfc, y=Family, yend=Family, color = association)) +
  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +ylab(NULL)# +

lollipop.family

## Genus ----
res.genus.sig.fix = res.genus.df.fix[res.genus.df.fix$q_value < 0.05,]  # separate out significant results

# Filter results to only include more prevalent taxa (otherwise lollipop plot is too busy)
sum_abundance = as.data.frame(rowSums(table))   #Create dataframe with OTU counts
res.genus.sig.fix$sum_abundance = match(res.genus.sig.fix$taxon, rownames(sum_abundance[1])) #Add counts to significant results dataframe
res.genus.sig.fix = res.genus.sig.fix[order(res.genus.sig.fix$sum_abundance, decreasing = TRUE),] #Sort by counts column (decreasing)
res.genus.sig.fix.30 = res.genus.sig.fix[1:30,] #Subset first 30 rows (i.e. first 30 most prevalent significantly differentially abundant taxa)

ggplot(res.genus.sig.fix.30, aes(lfc, reorder(Genus,sum_abundance), color = association ))+
  geom_point() +
  scale_y_discrete(position = "right") +
  theme_bw(base_size=15) +
  geom_segment(aes(x=0, xend=lfc, y=Genus, yend=Genus, color = association),size=1) +
  geom_vline(xintercept = 0, size=0.3) +xlab("log2FoldChange") +ylab(NULL)# +

ggsave("figures/ancom/lollipop-genus-30.svg", dpi=300)

# Volcano plot ----
## Phylum ----
res.phylum.df$association = ifelse(res.phylum.df$q_value < 0.05 & res.phylum.df$lfc > 0, "Lactating", "Not significant")
res.phylum.df$association = ifelse(res.phylum.df$q_value < 0.05 & res.phylum.df$lfc < 0, "Non-lactating",res.phylum.df$association)

ggplot(res.phylum.df, aes(x = as.numeric(lfc), y = -log10(as.numeric(q_value)), color = association)) +
  geom_point(size=3) +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_text_repel(data=res.phylum.sig, aes(label=Phylum), size=5, show.legend = FALSE) +
  theme_bw(base_size=20) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_manual(values=c("brown1", "darkcyan", "gray")) +
  theme(aspect.ratio=1)
  
ggsave("figures/ancom/volcano-phylum.svg", dpi=300)

## Family ----
res.family.df$association = ifelse(res.family.df$q_value < 0.05 & res.family.df$lfc > 0, "Lactating", "Not significant")
res.family.df$association = ifelse(res.family.df$q_value < 0.05 & res.family.df$lfc < 0, "Non-lactating",res.family.df$association)

ggplot(res.family.df, aes(x = as.numeric(lfc), y = -log10(as.numeric(q_value)), color = association)) +
  geom_point() +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot (Family)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  theme_bw()

## Genus ----
res.genus.df$association = ifelse(res.genus.df$q_value < 0.05 & res.genus.df$lfc > 0, "Lactating", "Not significant")
res.genus.df$association = ifelse(res.genus.df$q_value < 0.05 & res.genus.df$lfc < 0, "Non-lactating",res.genus.df$association)

ggplot(res.genus.df, aes(x = as.numeric(lfc), y = -log10(as.numeric(q_value)), color = association)) +
  geom_point() +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot (Genus)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_text_repel(data=res.genus.sig, aes(label=Genus),  max.overlaps = 18) +
  theme_bw()
