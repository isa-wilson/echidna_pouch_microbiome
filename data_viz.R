# Echidna pouch microbiome data visualisation
# Isabella Wilson 2024-12-20

library(tidyverse)
library(qiime2R)
library(phyloseq)
library(fantaxtic)
library(microViz)
library(microbiome)
library(svglite)

# Taxa bar plots ----
  
  ## Prepare data ----  
  
  metadata <- read_tsv("sample-metadata-pouch.tsv")
  metadata.df <- as.data.frame(metadata) # Convert to dataframe
  table <- read_qza("filtered-table-pouch.qza")$data
  taxonomy<-read_qza("taxonomy.qza")$data %>% parse_taxonomy()
  taxonomy.matrix <- as.matrix(taxonomy) # Convert to matrix
  
  # Create phyloseq object
  OTU = otu_table(table, taxa_are_rows=TRUE)
  TAX = tax_table(taxonomy.matrix)
  sampledata = sample_data(metadata.df)
  rownames(sampledata) <- metadata.df$`sample-id`
  
  physeq = phyloseq(OTU,TAX,sampledata)
  
  ## Generate top 15 genera bar plot  ----
    # Rename "NA"/"uncultured" taxa with last identified rank
    physeq.fix <- tax_fix(physeq, unknowns = c("Ambiguous_taxa", "uncultured", "NA", "Unknown_Family"))
  
    # Plot taxa                   
    physeq.fix %>%
      comp_barplot("Genus", n_taxa = 15, merge_other = FALSE, bar_outline_colour = "#f9f9f9", label = "CaptivityStatus", sample_order = c("2P","2B","7B","7PM","MSSBE114","MSSBE115","MSSBE116","MSSBE117","MSSBE118","MSSBE119","MSSBE113","NBSBE1","NBSBE2","NBSBE3","IWSBE5","1B","3B","5B","6B","8B","9B","10B")) +
      facet_grid(~factor(BreedingStatus, levels = c("Non-breeding", "Breeding", "Lactating")), scales = "free", space = "free") +
      theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5))
  
    # Export
    ggsave(filename = "figures/genus-plot.svg", dpi = 300)
    
  ## Generate nested plot ----
    # Obtain most abundant phyla and most abundant genera within those phyla
    top_nested <- nested_top_taxa(physeq,
                                  top_tax_level = "Phylum",
                                  nested_tax_level = "Genus",
                                  n_top_taxa = 4, 
                                  n_nested_taxa = 4)
    
    # Plot the relative abundances at two levels
    basic.plot <- plot_nested_bar(top_nested$ps_obj,
                                  top_level = "Phylum",
                                  nested_level = "Genus",
                                  nested_merged_label = "NA and other <tax>",
                                  legend_title = "Phylum and Genus",
                                  palette = c(Actinobacteriota = "blue", 
                                              Bacteroidota = "yellow", 
                                              Firmicutes = "pink", 
                                              Proteobacteria = "green"),
                                  sample_order = c("2P","2B","7B","7PM","MSSBE114","MSSBE115","MSSBE116","MSSBE117","MSSBE118","MSSBE119","MSSBE113","NBSBE1","NBSBE2","NBSBE3","IWSBE5","1B","3B","5B","6B","8B","9B","10B"))
    
    # Add faceting
    basic.plot +
      facet_grid(~factor(BreedingStatus, levels = c("Non-breeding", "Breeding", "Lactating")), scales = "free", space = "free") +
      guides(fill = guide_legend(ncol = 1)) +
      theme(axis.text.x = element_text(angle = 90)) + 
      theme(strip.background = element_rect(fill="lightgray")) +
      theme(strip.text = element_text(colour = 'black', size = 10, face = "bold")) +
      theme(legend.key.size = unit(15, "pt")) 
  
    # Export
    ggsave(filename = "figures/nested-plot.svg", dpi = 300)
  
    
# PCOA plots ----
    
  ## Captive vs Wild ----
    
    # Read in data
    taxonomy.PCOA <- read_qza("taxonomy.qza")$data %>% parse_taxonomy()
    taxonomy.PCOA <- as.matrix(taxonomy.PCOA) # convert to matrix
    rooted.tree <- read_tree("rooted-tree.nwk")
    table.captivity <- read_qza("filtered-table-pouch-non-breeding.qza")$data
    metadata.captivity <- read_tsv("sample-metadata-pouch-non-breeding.tsv")
    metadata.captivity <- as.data.frame(metadata.captivity) # convert to dataframe
    
    # Create physeq object
    OTU.captivity = otu_table(table.captivity, taxa_are_rows=TRUE)
    TAX.captivity = tax_table(taxonomy.PCOA)
    sampledata.captivity = sample_data(metadata.captivity)
    rownames(sampledata.captivity) <- metadata.captivity$`sample-id`
    
    physeq.captivity = phyloseq(OTU.captivity,TAX.captivity,sampledata.captivity,rooted.tree)
    
  ### Unweighted ----
      ordu.captivity = ordinate(physeq.captivity, "PCoA", "unifrac", weighted=FALSE)
      PCOA.captivity.unweighted.1.2 = plot_ordination(physeq.captivity, ordu.captivity, color="CaptivityStatus", axes = c(1,2))
      PCOA.captivity.unweighted.1.3 = plot_ordination(physeq.captivity, ordu.captivity, color="CaptivityStatus", axes = c(1,3))
      
      # Axes 1 and 2
      PCOA.captivity.unweighted.1.2 +
        scale_colour_manual(values = c("Captive" = "#f27304", "Wild" = "#008000")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )

      ggsave(filename = "figures/PCOA-captivity-unweighted-axes-1-2.svg", dpi = 300)
      
      # Axes 1 and 3
      PCOA.captivity.unweighted.1.3 +
        scale_colour_manual(values = c("Captive" = "#f27304", "Wild" = "#008000")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-captivity-unweighted-axes-1-3.svg", dpi = 300)

  ### Weighted ----
      ordu.captivity = ordinate(physeq.captivity, "PCoA", "unifrac", weighted=TRUE)
      PCOA.captivity.weighted.1.2 = plot_ordination(physeq.captivity, ordu.captivity, color="CaptivityStatus", axes = c(1,2))
      PCOA.captivity.weighted.1.3 = plot_ordination(physeq.captivity, ordu.captivity, color="CaptivityStatus", axes = c(1,3))
      
      # Axes 1 and 2
      PCOA.captivity.weighted.1.2 +
        scale_colour_manual(values = c("Captive" = "#f27304", "Wild" = "#008000")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-captivity-weighted-axes-1-2.svg", dpi = 300)
      
      # Axes 1 and 3
      PCOA.captivity.weighted.1.3 +
        scale_colour_manual(values = c("Captive" = "#f27304", "Wild" = "#008000")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-captivity-weighted-axes-1-3.svg", dpi = 300)
    
  ## Breeding vs non-breeding ----
        
    # Read in data
    table.breeding <- read_qza("filtered-table-pouch-non-lactating-wild.qza")$data
    metadata.breeding <- read_tsv("sample-metadata-pouch-non-lactating-wild.tsv")
    metadata.breeding <- as.data.frame(metadata.breeding) # convert to dataframe
    
    # Create physeq object
    OTU.breeding = otu_table(table.breeding, taxa_are_rows=TRUE)
    TAX.breeding = tax_table(taxonomy.PCOA)
    sampledata.breeding = sample_data(metadata.breeding)
    rownames(sampledata.breeding) <- metadata.breeding$`sample-id`
    
    physeq.breeding = phyloseq(OTU.breeding,TAX.breeding,sampledata.breeding,rooted.tree)
    
  ### Unweighted ----
      ordu.breeding = ordinate(physeq.breeding, "PCoA", "unifrac", weighted=FALSE)
      PCOA.breeding.unweighted.1.2 = plot_ordination(physeq.breeding, ordu.breeding, color="BreedingStatus", axes = c(1,2))
      PCOA.breeding.unweighted.1.3 = plot_ordination(physeq.breeding, ordu.breeding, color="BreedingStatus", axes = c(1,3))
      
      # Axes 1 and 2
      PCOA.breeding.unweighted.1.2 +
        scale_colour_manual(values = c("Breeding" = "purple", "Non-breeding" = "red")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-breeding-unweighted-axes-1-2.svg", dpi = 300)
      
      # Axes 1 and 3
      PCOA.breeding.unweighted.1.3 +
        scale_colour_manual(values = c("Breeding" = "purple", "Non-breeding" = "red")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-breeding-unweighted-axes-1-3.svg", dpi = 300)
      
  ### Weighted ----
      ordu.breeding = ordinate(physeq.breeding, "PCoA", "unifrac", weighted=TRUE)
      PCOA.breeding.weighted.1.2 = plot_ordination(physeq.breeding, ordu.breeding, color="BreedingStatus", axes = c(1,2))
      PCOA.breeding.weighted.1.3 = plot_ordination(physeq.breeding, ordu.breeding, color="BreedingStatus", axes = c(1,3))
      
      # Axes 1 and 2
      PCOA.breeding.weighted.1.2 +
        scale_colour_manual(values = c("Breeding" = "purple", "Non-breeding" = "red")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-breeding-weighted-axes-1-2.svg", dpi = 300)
      
      # Axes 1 and 3
      PCOA.breeding.weighted.1.3 +
        scale_colour_manual(values = c("Breeding" = "purple", "Non-breeding" = "red")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-breeding-weighted-axes-1-3.svg", dpi = 300)      

  ## Lactating vs non-lactating ----
        
    # Read in data
    table.lact <- read_qza("filtered-table-pouch.qza")$data
    metadata.lact <- read_tsv("sample-metadata-pouch.tsv")
    metadata.lact <- as.data.frame(metadata.lact) # convert to dataframe
    
    # Create physeq object
    OTU.lact = otu_table(table.lact, taxa_are_rows=TRUE)
    TAX.lact = tax_table(taxonomy.PCOA)
    sampledata.lact = sample_data(metadata.lact)
    rownames(sampledata.lact) <- metadata.lact$`sample-id`
    
    physeq.lact = phyloseq(OTU.lact,TAX.lact,sampledata.lact,rooted.tree)
    
  ### Unweighted ----
      ordu.lact = ordinate(physeq.lact, "PCoA", "unifrac", weighted=FALSE)
      PCOA.lact.unweighted.1.2 = plot_ordination(physeq.lact, ordu.lact, color="SampleType", axes = c(1,2))
      PCOA.lact.unweighted.1.3 = plot_ordination(physeq.lact, ordu.lact, color="SampleType", axes = c(1,3))
      
      # Axes 1 and 2
      PCOA.lact.unweighted.1.2 +
        scale_colour_manual(values = c("Lactating pouch" = "pink", "Non-lactating pouch" = "blue")) +
        stat_ellipse() +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-lactating-unweighted-axes-1-2.svg", dpi = 300)
      
      # Axes 1 and 3
      PCOA.lact.unweighted.1.3 +
        scale_colour_manual(values = c("Lactating pouch" = "pink", "Non-lactating pouch" = "blue")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-lactating-unweighted-axes-1-3.svg", dpi = 300)
      
  ### Weighted ----
      ordu.lact = ordinate(physeq.lact, "PCoA", "unifrac", weighted=TRUE)
      PCOA.lact.weighted.1.2 = plot_ordination(physeq.lact, ordu.lact, color="SampleType", axes = c(1,2))
      PCOA.lact.weighted.1.3 = plot_ordination(physeq.lact, ordu.lact, color="SampleType", axes = c(1,3))
      
      # Axes 1 and 2
      PCOA.lact.weighted.1.2 +
        scale_colour_manual(values = c("Lactating pouch" = "pink", "Non-lactating pouch" = "blue")) +
        stat_ellipse() +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-lactating-weighted-axes-1-2.svg", dpi = 300)
      
      # Axes 1 and 3
      PCOA.lact.weighted.1.3 +
        scale_colour_manual(values = c("Lactating pouch" = "pink", "Non-lactating pouch" = "blue")) +
        geom_point(size = 3.5) +
        theme_bw() +
        theme(panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black", linewidth = 0.75),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_text(size = 12),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)
              )
      
      ggsave(filename = "figures/PCOA-lactating-weighted-axes-1-3.svg", dpi = 300) 
    
# Comparing all sample types ----
  ## Prepare data ----
    table.all <- read_qza("freq-decontam-filtered-table.qza")$data
    metadata.all <- read_tsv("sample-metadata-all.tsv")
    metadata.all <- as.data.frame(metadata.all) # convert to dataframe
    
    # Create phyloseq object
    OTU.all = otu_table(table.all, taxa_are_rows=TRUE)
    sampledata.all = sample_data(metadata.all)
    rownames(sampledata.all) <- metadata.all$`sample-id`
    
    physeq.all = phyloseq(OTU.all,TAX,sampledata.all,rooted.tree)
  
  ## Plot PCOA ----
    ### Unweighted ----
    ordu.all = ordinate(physeq.all, "PCoA", "unifrac", weighted=FALSE)
    PCOA.all.unweighted.1.2 = plot_ordination(physeq.all, ordu.all, color="SampleType", axes = c(1,2))
    PCOA.all.unweighted.1.3 = plot_ordination(physeq.all, ordu.all, color="SampleType", axes = c(1,3))
    
    # Axes 1 and 2
    PCOA.all.unweighted.1.2 +
      scale_colour_manual(values = c("Lactating pouch" = "orchid1", "Non-lactating pouch" = "blue", "Cloacal" = "red", "Oral" = "deepskyblue","Environment" = "chartreuse3", "Negative control" = "darkgoldenrod1" )) +
      geom_point(size = 2.5) +
      theme_bw() +
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black", linewidth = 0.75),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_blank(),
            axis.title.y = element_text(size = 12),
            legend.position = "right",
            legend.text = element_text(size = 12)
            )
    
    ggsave(filename = "figures/PCOA-all-unweighted-axes-1-2.svg", dpi = 300)
    
    # Axes 1 and 3
    PCOA.all.unweighted.1.3 +
      scale_colour_manual(values = c("Lactating pouch" = "orchid1", "Non-lactating pouch" = "blue", "Cloacal" = "red", "Oral" = "deepskyblue","Environment" = "chartreuse3", "Negative control" = "darkgoldenrod1" )) +
      geom_point(size = 2.5) +
      theme_bw() +
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black", linewidth = 0.75),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_blank(),
            axis.title.y = element_text(size = 12),
            legend.position = "right",
            legend.text = element_text(size = 12)
            )
    ggsave(filename = "figures/PCOA-all-unweighted-axes-1-3.svg", dpi = 300)
    
    ### Weighted ----
    ordu.all = ordinate(physeq.all, "PCoA", "unifrac", weighted=TRUE)
    PCOA.all.weighted.1.2 = plot_ordination(physeq.all, ordu.all, color="SampleType", axes = c(1,2))
    PCOA.all.weighted.1.3 = plot_ordination(physeq.all, ordu.all, color="SampleType", axes = c(1,3))
    
    # Axes 1 and 2
    PCOA.all.weighted.1.2 +
      scale_colour_manual(values = c("Lactating pouch" = "orchid1", "Non-lactating pouch" = "blue", "Cloacal" = "red", "Oral" = "deepskyblue","Environment" = "chartreuse3", "Negative control" = "darkgoldenrod1" )) +
      geom_point(size = 2.5) +
      theme_bw() +
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black", linewidth = 0.75),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_blank(),
            axis.title.y = element_text(size = 12),
            legend.position = "right",
            legend.text = element_text(size = 12)
            )
    
    ggsave(filename = "figures/PCOA-all-weighted-axes-1-2.svg", dpi = 300)
    
    # Axes 1 and 3
    PCOA.all.weighted.1.3 +
      scale_colour_manual(values = c("Lactating pouch" = "orchid1", "Non-lactating pouch" = "blue", "Cloacal" = "red", "Oral" = "deepskyblue","Environment" = "chartreuse3", "Negative control" = "darkgoldenrod1" )) +
      geom_point(size = 2.5) +
      theme_bw() +
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black", linewidth = 0.75),
            axis.ticks = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.text.y = element_blank(),
            axis.title.y = element_text(size = 12),
            legend.position = "right",
            legend.text = element_text(size = 12)
            )
    
    ggsave(filename = "figures/PCOA-all-weighted-axes-1-3.svg", dpi = 300)
