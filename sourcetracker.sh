# Identifying the origins of echidna pouch microbes using Sourcetracker2

# Isabella Wilson 2024-12-20
# Sourcetracker2 version 2.0.1-dev

source activate st2

sourcetracker2 gibbs \
	-i feature-table.biom \
	-m sourcetracker-metadata.txt \
	--source_sink_column SourceSink \
    --source_column_value source \
    --sink_column_value sink \
    --source_category_column Env \
    --sink_rarefaction_depth 6800 \
    --per_sink_feature_assignments \
    -o sourcetracker_results_OTUs
