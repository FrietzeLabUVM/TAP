scripts=~/lab_shared/scripts/vacc_rnaseq_pipeline
dock_test_dir=$scripts/testing_docker
cd $dock_test_dir
test_dir=~/lab_shared/scripts/vacc_rnaseq_pipeline/testing
#args2="-noSub"
args3="--docker tap"
args="-i $test_dir/test_data -ref $(readlink -f ~/lab_shared/indexes/DM6)"

#no config file
bash $scripts/submit_rnaseq_pipeline.sh -o $dock_test_dir/test_alignment.no_config $args $args2 ${args3}

exit 0
#config that just runs sample 2
bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.basic.csv -o $dock_test_dir/test_alignment.basic $args $args2 ${args3}

#config that renames samples
bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.rename.csv -o $dock_test_dir/test_alignment.rename $args $args2 ${args3}

#config that pools sample 1 and 2, leaves 3 alone, then pools all using & and " "
bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.pool.csv -o $dock_test_dir/test_alignment.pool $args $args2 ${args3}

#config that controls all params, config is only argument
bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.params.csv $args2 ${args3}

#config for SE mode, one pooled too
bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.SE.csv $args2 ${args3}

#config for running rDNA only
bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.rDNA_only.csv $args2 ${args3}

#config for running ChIP-seq
bash $scripts/submit_chipseq_pipeline.sh -c $test_dir/test_data/chip_test_hg38_config.params.csv ${args3}

#config for running ChIP-seq using pooling
bash $scripts/submit_chipseq_pipeline.sh -c $test_dir/test_data/chip_test_hg38_config.pool_params.csv ${args3}
