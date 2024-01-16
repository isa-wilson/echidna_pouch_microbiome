### QIIME2 WORKFLOW ###

# Isabella Wilson 02/05/2023
# Qiime2 version 2020.2

##########################################

#Activate qiime2
conda activate qiime2-2020.2


# Import manifest file

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path /path/to/manifest.tsv \
  --output-path paired-end-demux.qza

	#view demux sequences
	qiime demux summarize \
	  --i-data paired-end-demux.qza \
	  --o-visualization demux_seqs.qzv
	  
	  
# Join paired reads

qiime vsearch join-pairs \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --o-joined-sequences demux-joined.qza
  
  #view paired reads
	qiime demux summarize \
	  --i-data demux-joined.qza \
	  --o-visualization demux-joined.qzv
	  
	  
# Quality filter joined reads
	  
qiime quality-filter q-score-joined \
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
	  
	# Run decontam - see decontam_script.R file
	
	# Create table excluding contaminants
	qiime feature-table filter-features \
		--i-table freq-filtered-table.qza \
		--p-exclude-ids \
		--m-metadata-file contaminant_features.txt \
		--o-filtered-table freq-decontam-filtered-table.qza
		
	# Visualise table
	qiime feature-table summarize \
	  --i-table freq-decontam-filtered-table.qza \
	  --o-visualization freq-decontam-filtered-table.qzv
 	 
# Generate tree for phylogenetic diversity analysis

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
  
# Alpha rarefaction (max depth based on minimum frequency from table)
	# Note - sample 9B had an unexpectedly low sampling depth (lower than negative controls)
	
qiime diversity alpha-rarefaction \
  --i-table freq-decontam-filtered-table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 6300 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
  
# Core metrics alpha and beta diversity (sampling depth based on alpha rarefaction curve)

# Analysing all samples

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table freq-decontam-filtered-table.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results
  
  # Generate alpha diversity plots
  
	# Faith's PD
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
	  --m-metadata-file sample-metadata.tsv \
	  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
	  
	# Pielou's evenness
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
	  --m-metadata-file sample-metadata.tsv \
	  --o-visualization core-metrics-results/evenness-group-significance.qzv
	  
  	# Shannon diversity
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
	  --m-metadata-file sample-metadata.tsv \
	  --o-visualization core-metrics-results/shannon-group-significance.qzv
	  
	# Observed OTUs
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
	  --m-metadata-file sample-metadata.tsv \
	  --o-visualization core-metrics-results/observed_otus.qzv
	  
  # Generate beta diversity stats
	
	# For sample type  
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata.tsv \
	  --m-metadata-column SampleType \
	  --o-visualization core-metrics-results/unweighted-unifrac-sample-type-significance.qzv \
	  --p-pairwise

	# For reproductive status
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata.tsv \
	  --m-metadata-column BreedingStatus \
	  --o-visualization core-metrics-results/unweighted-unifrac-breeding-status-significance.qzv \
	  --p-pairwise

	# For captivity status
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata.tsv \
	  --m-metadata-column CaptivityStatus \
	  --o-visualization core-metrics-results/unweighted-unifrac-captivity-significance.qzv \
	  --p-pairwise


# Analysing pouch samples only

qiime feature-table filter-samples \
	  --i-table freq-decontam-filtered-table.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --o-filtered-table filtered-table-pouch.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table filtered-table-pouch.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file sample-metadata-pouch.tsv \
  --output-dir core-metrics-results-pouch
  
  # Generate alpha diversity plots
  
	# Faith's PD
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results-pouch/faith_pd_vector.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --o-visualization core-metrics-results-pouch/faith-pd-group-significance.qzv
	  
	# Pielou's evenness
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results-pouch/evenness_vector.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --o-visualization core-metrics-results-pouch/evenness-group-significance.qzv
	  
  	# Shannon diversity
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results-pouch/shannon_vector.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --o-visualization core-metrics-results-pouch/shannon-group-significance.qzv
	  
	# Observed OTUs
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results-pouch/observed_features_vector.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --o-visualization core-metrics-results-pouch/observed_otus.qzv
	  
  # Generate beta diversity stats
	
	# For sample type  
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results-pouch/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --m-metadata-column SampleType \
	  --o-visualization core-metrics-results-pouch/unweighted-unifrac-sample-type-significance.qzv \
	  --p-pairwise

	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results-pouch/weighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --m-metadata-column SampleType \
	  --o-visualization core-metrics-results-pouch/weighted-unifrac-sample-type-significance.qzv \
	  --p-pairwise
	  
	# For reproductive status
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results-pouch/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --m-metadata-column BreedingStatus \
	  --o-visualization core-metrics-results-pouch/unweighted-unifrac-breeding-status-significance.qzv \
	  --p-pairwise

	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results-pouch/weighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --m-metadata-column BreedingStatus \
	  --o-visualization core-metrics-results-pouch/weighted-unifrac-breeding-status-significance.qzv \
	  --p-pairwise
	  
	# For captivity status
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results-pouch/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --m-metadata-column CaptivityStatus \
	  --o-visualization core-metrics-results-pouch/unweighted-unifrac-captivity-significance.qzv \
	  --p-pairwise

	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results-pouch/weighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata-pouch.tsv \
	  --m-metadata-column CaptivityStatus \
	  --o-visualization core-metrics-results-pouch/weighted-unifrac-captivity-significance.qzv \
	  --p-pairwise

# Analysing non-lactating samples only (for captivity analysis)

qiime feature-table filter-samples \
	  --i-table freq-decontam-filtered-table.qza \
	  --m-metadata-file sample-metadata-pouch-nonlactating.tsv \
	  --o-filtered-table filtered-table-pouch-nonlactating.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table filtered-table-pouch-nonlactating.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file sample-metadata-pouch-nonlactating.tsv \
  --output-dir core-metrics-results-pouch-nonlactating
  
  # Generate alpha diversity plots
  
	# Faith's PD
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results-pouch-nonlactating/faith_pd_vector.qza \
	  --m-metadata-file sample-metadata-pouch-nonlactating.tsv \
	  --o-visualization core-metrics-results-pouch-nonlactating/faith-pd-group-significance.qzv
	  
	# Pielou's evenness
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results-pouch-nonlactating/evenness_vector.qza \
	  --m-metadata-file sample-metadata-pouch-nonlactating.tsv \
	  --o-visualization core-metrics-results-pouch-nonlactating/evenness-group-significance.qzv
	  
  	# Shannon diversity
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results-pouch-nonlactating/shannon_vector.qza \
	  --m-metadata-file sample-metadata-pouch-nonlactating.tsv \
	  --o-visualization core-metrics-results-pouch-nonlactating/shannon-group-significance.qzv
	  
	# Observed OTUs
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results-pouch-nonlactating/observed_features_vector.qza \
	  --m-metadata-file sample-metadata-pouch-nonlactating.tsv \
	  --o-visualization core-metrics-results-pouch-nonlactating/observed_otus.qzv
	  	 
  # Generate beta diversity stats	 
  
	# Unweighted Unifrac  
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results-pouch-nonlactating/unweighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata-pouch-nonlactating.tsv \
	  --m-metadata-column CaptivityStatus \
	  --o-visualization core-metrics-results-pouch-nonlactating/unweighted-unifrac-captivity-significance.qzv \
	  --p-pairwise

	# Weighted Unifrac  
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results-pouch-nonlactating/weighted_unifrac_distance_matrix.qza \
	  --m-metadata-file sample-metadata-pouch-nonlactating.tsv \
	  --m-metadata-column CaptivityStatus \
	  --o-visualization core-metrics-results-pouch-nonlactating/weighted-unifrac-captivity-significance.qzv \
	  --p-pairwise
	  
	  
# Classify taxa

qiime feature-classifier classify-sklearn \
  --i-classifier SILVA-v138-515f-806r-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Generate barplots  

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