#!/bin/bash

/home/emoeller/anaconda3/envs/scBonita/bin/python pipeline.py \
    --data_file "../network_simulation/boolean_network_conditions/network_simulation.csv" \
    --dataset_name "boolean" \
    --full_pipeline 1 \
    --pathway_list "../network_simulation/boolean_network_conditions/random_boolean_network.graphml" \
    --max_nodes 20000 \
    --separator ',' \
    --get_kegg_pathways False \
    --network "../network_simulation/boolean_network_conditions/random_boolean_network.graphml" \
    --organism "hsa" \
    --conda_env "scBonita" \
    --python_version "/home/emoeller/anaconda3/envs/scBonita/bin/python" \
    --binarize_threshold 0.1 \
    --display_title "True" \
    --run_all_networks "True" \
    --run_attractor_analysis "True"
