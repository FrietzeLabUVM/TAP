
#no config file
#bash submit_rnaseq_pipeline.sh -i test_data -o test_alignment.no_config -ref references/DM6

#config that just runs sample 2
#bash submit_rnaseq_pipeline.sh -c test_data/test_dm6_config.basic.csv -i test_data -o test_alignment.basic -ref references/DM6

#config that pools sample 1 and 2, leaves 3 alone
bash submit_rnaseq_pipeline.sh -c test_data/test_dm6_config.pool.csv -i test_data -o test_alignment.pool -ref references/DM6

