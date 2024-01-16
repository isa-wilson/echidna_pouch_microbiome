### SOURCETRACKER ###

# Isabella Wilson 02/05/2023
# sourcetracker2, version 2.0.1.dev0

##########################################

# Before starting, create metadata file
# 	1) Create SinkSource column; label all sources as "source", all sinks as "sink"
# 	2) Rename sample-type column to Env
	
# Export Qiime OTU table as BIOM table
conda activate qiime2-2020.2

qiime tools export \
 --input-path freq-decontam-filtered-table.qza \
 --output-path exported-feature-table

#Run Sourcetracker2
source activate st2

sourcetracker2 gibbs \
    -i feature-table.biom \
    -m sourcetracker-metadata.txt \
    --source_sink_column SourceSink \
    --source_column_value source \
    --sink_column_value sink \
    --source_category_column Env \
    --sink_rarefaction_depth 1000 \
    -o sourcetracker_results
    
#Track contribution of each feature to the sources    
sourcetracker2 \
	-i feature-table.biom \
	-m sourcetracker-metadata.txt \
	--source_sink_column SourceSink \
    --source_column_value source \
    --sink_column_value sink \
    --source_category_column Env \
    --sink_rarefaction_depth 1000 \
    --per_sink_feature_assignments \
    -o sourcetracker_results_OTUs