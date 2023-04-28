set -e
dock_test_dir=$(dirname "$0")
dock_test_dir=$(readlink -f $dock_test_dir)
scripts=$(dirname $dock_test_dir)
#dock_test_dir=$scripts/testing_docker
#cd $dock_test_dir
test_dir=${scripts}/testing
#args2="-noSub"
args3="--docker jrboyd/tap -noSub"
args3="--singularity tap_latest.sif -noSub"
args="-i $test_dir/test_data -ref $(readlink -f ~/lab_shared/indexes/DM6)"

#no config file
#bash $scripts/submit_rnaseq_pipeline.sh -o $dock_test_dir/test_alignment.no_config -SE $args $args2 ${args3}


#config that just runs sample 2
#bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.basic.csv -o $dock_test_dir/test_alignment.basic $args $args2 ${args3}

#config that renames samples
#bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.rename.csv -o $dock_test_dir/test_alignment.rename $args $args2 ${args3}

#config that pools sample 1 and 2, leaves 3 alone, then pools all using & and " "
#bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.pool.csv -o $dock_test_dir/test_alignment.pool $args $args2 ${args3}

#config that controls all params, config is only argument
#bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.params.csv $args2 ${args3}

#config for SE mode, one pooled too
#bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.SE.csv $args2 ${args3}

#config for running rDNA only
#bash $scripts/submit_rnaseq_pipeline.sh -c $test_dir/test_data/test_dm6_config.rDNA_only.csv $args2 ${args3}

#config for running ChIP-seq
bash $scripts/submit_chipseq_pipeline.sh -c $test_dir/test_data/chip_test_hg38_config.params.csv ${args3}

#config for running ChIP-seq using pooling
bash $scripts/submit_chipseq_pipeline.sh -c $test_dir/test_data/chip_test_hg38_config.pool_params.csv ${args3}

