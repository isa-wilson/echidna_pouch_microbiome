### QIIME2 WORKFLOW ###

# Isabella Wilson 2024-12-20
# Qiime2 version 2024.10

##########################################

#Activate qiime2
conda activate qiime2-amplicon-2024.10

# Import manifest file

    qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-format PairedEndFastqManifestPhred33V2 \
    --input-path /Users/peterconway/bella_working/qiime4/manifest_monotreme.tsv \
    --output-path paired-end-demux.qza

	#view demux sequences
	qiime demux summarize \
	--i-data paired-end-demux.qza \
	--o-visualization demux_seqs.qzv
	  
	  
# Join paired reads

    qiime vsearch merge-pairs \
    --i-demultiplexed-seqs paired-end-demux.qza \
    --o-merged-sequences demux-joined.qza \
    --o-unmerged-sequences demux-unmerged.qza

    #view paired reads
    qiime demux summarize \
    --i-data demux-joined.qza \
    --o-visualization demux-joined.qzv
	  
	  
# Quality filter joined reads
	  
    qiime quality-filter q-score \
    --i-demux demux-joined.qza \
    --o-filtered-sequences demux-joined-filtered.qza \
    --o-filter-stats demux-joined-filter-stats.qza
    
    #view quality filter stats
	qiime metadata tabulate \
 	--m-input-file demux-joined-filter-stats.qza \
 	--o-visualization demux-joined-filter-stats.qzv
  
# Use deblur and trim based on visualising joined reads
  
    qiime deblur denoise-16S \
    --i-demultiplexed-seqs demux-joined-filtered.qza \
    --p-trim-length 252 \
    --p-sample-stats \
    --o-representative-sequences rep-seqs.qza \
    --o-table table.qza \
    --o-stats deblur-stats.qza

 	# Remove singletons
 	qiime feature-table filter-features \
	--i-table table.qza \
	--p-min-frequency 2 \
	--o-filtered-table freq-filtered-table.qza
	  
# Generate tree for phylogenetic diversity analysis

    qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs.qza \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza

    qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads SILVA/silva-138-99-seqs-515-806.qza \
    --i-reference-taxonomy SILVA/silva-138-99-tax-515-806.qza \
    --o-classifier classifier.qza

# Classify taxa
    qiime feature-classifier classify-sklearn \
    --i-classifier classifier.qza \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza

    qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv
	
### RUN DECONTAM - see decontam_script.R file ###

# Exclude contaminants from feature table
	
    # Remove decontam-identified contaminants
	qiime feature-table filter-features \
	--i-table freq-filtered-table.qza \
    --p-exclude-ids \
	--m-metadata-file decontam/contaminant.ASVs.txt \
	--o-filtered-table freq-decontam-table.qza

    # Remove mitochondria and chloroplast
	qiime taxa filter-table \
	--i-table freq-decontam-table.qza \
    --i-taxonomy taxonomy.qza \
	--p-exclude mitochondria,chloroplast \
	--o-filtered-table freq-decontam-filtered-table.qza

	# Visualise table
	qiime feature-table summarize \
	--i-table freq-decontam-filtered-table.qza \
	--o-visualization freq-decontam-filtered-table.qzv

# Alpha rarefaction (max depth based on median frequency from table)
	# Note - sample 9B had an unexpectedly low sampling depth (lower than negative controls)
    
    qiime diversity alpha-rarefaction \
    --i-table freq-decontam-filtered-table.qza \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth 35000 \
    --m-metadata-file sample-metadata.tsv \
    --o-visualization alpha-rarefaction.qzv
  
# Perform alpha and beta diversity analysis (sampling depth selected using alpha rarefaction curve)

    # All samples
        qiime diversity core-metrics-phylogenetic \
        --i-phylogeny rooted-tree.qza \
        --i-table freq-decontam-filtered-table.qza \
        --p-sampling-depth 6800 \
        --m-metadata-file sample-metadata-all.tsv \
        --output-dir core-metrics-results
            
        # Generate beta diversity stats
            # Unweighted
            qiime diversity beta-group-significance \
            --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
            --m-metadata-file sample-metadata-all.tsv \
            --m-metadata-column SampleType \
            --o-visualization core-metrics-results/unweighted-unifrac-sample-type-significance.qzv \
            --p-pairwise
            
            # Weighted
            qiime diversity beta-group-significance \
            --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
            --m-metadata-file sample-metadata-all.tsv \
            --m-metadata-column SampleType \
            --o-visualization core-metrics-results/weighted-unifrac-sample-type-significance.qzv \
            --p-pairwise

    # Captive vs wild (non-breeding only)

        qiime feature-table filter-samples \
        --i-table freq-decontam-filtered-table.qza \
        --m-metadata-file sample-metadata-pouch-non-breeding.tsv \
        --o-filtered-table filtered-table-pouch-non-breeding.qza

        qiime diversity core-metrics-phylogenetic \
        --i-phylogeny rooted-tree.qza \
        --i-table filtered-table-pouch-non-breeding.qza \
        --p-sampling-depth 6800 \
        --m-metadata-file sample-metadata-pouch-non-breeding.tsv \
        --output-dir core-metrics-results-captivity
        
        # Generate alpha diversity plots
            # Faith's PD
            qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results-captivity/faith_pd_vector.qza \
            --m-metadata-file sample-metadata-pouch-non-breeding.tsv \
            --o-visualization core-metrics-results-captivity/faith-pd-group-significance.qzv
            
            # Pielou's evenness
            qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results-captivity/evenness_vector.qza \
            --m-metadata-file sample-metadata-pouch-non-breeding.tsv \
            --o-visualization core-metrics-results-captivity/evenness-group-significance.qzv
            
            # Shannon diversity
            qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results-captivity/shannon_vector.qza \
            --m-metadata-file sample-metadata-pouch-non-breeding.tsv \
            --o-visualization core-metrics-results-captivity/shannon-group-significance.qzv
            
            # Observed OTUs
            qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results-captivity/observed_features_vector.qza \
            --m-metadata-file sample-metadata-pouch-non-breeding.tsv \
            --o-visualization core-metrics-results-captivity/observed_otus.qzv
            
        # Generate beta diversity stats  
            
            qiime diversity beta-group-significance \
            --i-distance-matrix core-metrics-results-captivity/unweighted_unifrac_distance_matrix.qza \
            --m-metadata-file sample-metadata-pouch-non-breeding.tsv \
            --m-metadata-column CaptivityStatus \
            --o-visualization core-metrics-results-captivity/unweighted-unifrac-captivity-significance.qzv \
            --p-pairwise

            qiime diversity beta-group-significance \
            --i-distance-matrix core-metrics-results-captivity/weighted_unifrac_distance_matrix.qza \
            --m-metadata-file sample-metadata-pouch-non-breeding.tsv \
            --m-metadata-column CaptivityStatus \
            --o-visualization core-metrics-results-captivity/weighted-unifrac-captivity-significance.qzv \
            --p-pairwise

    # Breeding vs non-breeding (wild only)

        qiime feature-table filter-samples \
        --i-table freq-decontam-filtered-table.qza \
        --m-metadata-file sample-metadata-pouch-non-lactating-wild.tsv \
        --o-filtered-table filtered-table-pouch-non-lactating-wild.qza

        qiime diversity core-metrics-phylogenetic \
        --i-phylogeny rooted-tree.qza \
        --i-table filtered-table-pouch-non-lactating-wild.qza \
        --p-sampling-depth 6800 \
        --m-metadata-file sample-metadata-pouch-non-lactating-wild.tsv \
        --output-dir core-metrics-results-pouch-breeding
        
        # Generate alpha diversity plots
            # Faith's PD
            qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results-pouch-breeding/faith_pd_vector.qza \
            --m-metadata-file sample-metadata-pouch-non-lactating-wild.tsv \
            --o-visualization core-metrics-results-pouch-breeding/faith-pd-group-significance.qzv
            
            # Pielou's evenness
            qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results-pouch-breeding/evenness_vector.qza \
            --m-metadata-file sample-metadata-pouch-non-lactating-wild.tsv \
            --o-visualization core-metrics-results-pouch-breeding/evenness-group-significance.qzv
            
            # Shannon diversity
            qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results-pouch-breeding/shannon_vector.qza \
            --m-metadata-file sample-metadata-pouch-non-lactating-wild.tsv \
            --o-visualization core-metrics-results-pouch-breeding/shannon-group-significance.qzv
            
            # Observed OTUs
            qiime diversity alpha-group-significance \
            --i-alpha-diversity core-metrics-results-pouch-breeding/observed_features_vector.qza \
            --m-metadata-file sample-metadata-pouch-non-lactating-wild.tsv \
            --o-visualization core-metrics-results-pouch-breeding/observed_otus.qzv
                
        # Generate beta diversity stats	 
            # Unweighted Unifrac  
            qiime diversity beta-group-significance \
            --i-distance-matrix core-metrics-results-pouch-breeding/unweighted_unifrac_distance_matrix.qza \
            --m-metadata-file sample-metadata-pouch-non-lactating-wild.tsv \
            --m-metadata-column BreedingStatus \
            --o-visualization core-metrics-results-pouch-breeding/unweighted-unifrac-captivity-significance.qzv \
            --p-pairwise

            # Weighted Unifrac  
            qiime diversity beta-group-significance \
            --i-distance-matrix core-metrics-results-pouch-breeding/weighted_unifrac_distance_matrix.qza \
            --m-metadata-file sample-metadata-pouch-non-lactating-wild.tsv \
            --m-metadata-column BreedingStatus \
            --o-visualization core-metrics-results-pouch-breeding/weighted-unifrac-captivity-significance.qzv \
            --p-pairwise

    # Lactating vs non-lactating

            qiime feature-table filter-samples \
            --i-table freq-decontam-filtered-table.qza \
            --m-metadata-file sample-metadata-pouch.tsv \
            --o-filtered-table filtered-table-pouch.qza

            qiime diversity core-metrics-phylogenetic \
            --i-phylogeny rooted-tree.qza \
            --i-table filtered-table-pouch.qza \
            --p-sampling-depth 6800 \
            --m-metadata-file sample-metadata-pouch.tsv \
            --output-dir core-metrics-results-lactation
            
            # Generate alpha diversity plots
                # Faith's PD
                qiime diversity alpha-group-significance \
                --i-alpha-diversity core-metrics-results-lactation/faith_pd_vector.qza \
                --m-metadata-file sample-metadata-pouch.tsv \
                --o-visualization core-metrics-results-lactation/faith-pd-group-significance.qzv
                
                # Pielou's evenness
                qiime diversity alpha-group-significance \
                --i-alpha-diversity core-metrics-results-lactation/evenness_vector.qza \
                --m-metadata-file sample-metadata-pouch.tsv \
                --o-visualization core-metrics-results-lactation/evenness-group-significance.qzv
                
                # Shannon diversity
                qiime diversity alpha-group-significance \
                --i-alpha-diversity core-metrics-results-lactation/shannon_vector.qza \
                --m-metadata-file sample-metadata-pouch.tsv \
                --o-visualization core-metrics-results-lactation/shannon-group-significance.qzv
                
                # Observed OTUs
                qiime diversity alpha-group-significance \
                --i-alpha-diversity core-metrics-results-lactation/observed_features_vector.qza \
                --m-metadata-file sample-metadata-pouch.tsv \
                --o-visualization core-metrics-results-lactation/observed_otus.qzv
                    
            # Generate beta diversity stats	 
                # Unweighted Unifrac  
                qiime diversity beta-group-significance \
                --i-distance-matrix core-metrics-results-lactation/unweighted_unifrac_distance_matrix.qza \
                --m-metadata-file sample-metadata-pouch.tsv \
                --m-metadata-column SampleType \
                --o-visualization core-metrics-results-lactation/unweighted-unifrac-captivity-significance.qzv \
                --p-pairwise

                # Weighted Unifrac  
                qiime diversity beta-group-significance \
                --i-distance-matrix core-metrics-results-lactation/weighted_unifrac_distance_matrix.qza \
                --m-metadata-file sample-metadata-pouch.tsv \
                --m-metadata-column SampleType \
                --o-visualization core-metrics-results-lactation/weighted-unifrac-captivity-significance.qzv \
                --p-pairwise
# Generate taxonomy barplots  

	# All samples
        qiime taxa barplot \
        --i-table freq-decontam-filtered-table.qza \
        --i-taxonomy taxonomy.qza \
        --m-metadata-file sample-metadata.tsv \
        --o-visualization taxa-bar-plot-all.qzv

	# Pouch only
        qiime taxa barplot \
        --i-table filtered-table-pouch.qza \
        --i-taxonomy taxonomy.qza \
        --m-metadata-file sample-metadata-pouch.tsv \
        --o-visualization taxa-bar-plot-pouch.qzv
  
	# Grouping sample types for taxonomic analysis
        qiime feature-table group \
        --i-table filtered-table-pouch.qza \
        --p-axis sample \
        --m-metadata-file sample-metadata-pouch.tsv \
        --m-metadata-column BreedingStatus \
        --p-mode mean-ceiling \
        --o-grouped-table filtered-table-pouch-grouped.qza
        
        qiime taxa barplot \
        --i-table filtered-table-pouch-grouped.qza \
        --i-taxonomy taxonomy.qza \
        --m-metadata-file sample-metadata-pouch-grouped.tsv \
        --o-visualization taxa-bar-plot-pouch-grouped.qzv
	
	# Contaminants
        qiime taxa barplot \
        --i-table contaminant-filtered-table.qza \
        --i-taxonomy taxonomy.qza \
        --m-metadata-file sample-metadata.tsv \
        --o-visualization taxa-bar-plot-contaminant.qzv
	  
	# Negative controls (pre-decontam)
        qiime feature-table filter-samples \
        --i-table freq-filtered-table.qza \
        --m-metadata-file sample-metadata-neg.tsv \
        --o-filtered-table freq-filtered-table-neg.qza
        
        qiime taxa barplot \
        --i-table freq-filtered-table-neg.qza \
        --i-taxonomy taxonomy.qza \
        --m-metadata-file sample-metadata-neg.tsv \
        --o-visualization taxa-bar-plot-neg.qzv
	
	# Negative controls (post-decontam)
        qiime feature-table filter-samples \
        --i-table freq-decontam-filtered-table.qza \
        --m-metadata-file sample-metadata-neg.tsv \
        --o-filtered-table freq-decontam-filtered-table-neg.qza
        
        qiime taxa barplot \
        --i-table freq-decontam-filtered-table-neg.qza \
        --i-taxonomy taxonomy.qza \
        --m-metadata-file sample-metadata-neg.tsv \
        --o-visualization taxa-bar-plot-decontam-neg.qzv
