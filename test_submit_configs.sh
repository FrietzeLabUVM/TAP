
args2="-noSub"
args="-i test_data -ref references/DM6"

#no config file
#bash submit_rnaseq_pipeline.sh -o test_alignment.no_config $args $args2

#config that just runs sample 2
#bash submit_rnaseq_pipeline.sh -c test_data/test_dm6_config.basic.csv -o test_alignment.basic $args $args2

#config that renames samples
#bash submit_rnaseq_pipeline.sh -c test_data/test_dm6_config.rename.csv -o test_alignment.rename $args $args2

#config that pools sample 1 and 2, leaves 3 alone
#bash submit_rnaseq_pipeline.sh -c test_data/test_dm6_config.pool.csv -o test_alignment.pool $args $args2

#config that controls all params, config is only argument
bash submit_rnaseq_pipeline.sh -c test_data/test_dm6_config.params.csv $args2
