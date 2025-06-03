#!/bin/bash -l

# Base name (constant) for all runs
DATASET_NAME="simulated_data"

# List of gene‐counts to iterate over
GENE_LIST=(10 20 30 40 50 60 70 80 90 100)

for NUM_GENES in "${GENE_LIST[@]}"; do
    echo "========================================"
    echo "  Running pipeline for ${NUM_GENES} genes"
    echo "========================================"

    # 1) Generate the simulated network+data for this gene count
    python3 "/home/emoeller/github/scBONITA2/network_simulation/network_simulation.py" \
        --dataset_name "${DATASET_NAME}" \
        --num_genes "${NUM_GENES}"

    # 2) Run the main scBONITA pipeline
    python3 "/home/emoeller/github/scBONITA2/code/pipeline_class.py" \
        --data_file "input/${DATASET_NAME}_${NUM_GENES}.csv" \
        --dataset_name "${DATASET_NAME}_${NUM_GENES}" \
        --datafile_sep "," \
        --network_files "${DATASET_NAME}_${NUM_GENES}.graphml" \
        --binarize_threshold 0.01 \
        --get_kegg_pathways "False" \
        --minimum_overlap 1

    # 3) Compare scBONITA’s learned rules vs. the “truth” rules
    python3 "/home/emoeller/github/scBONITA2/network_simulation/check_rules_output.py" \
        --scbonita_rules "/home/emoeller/github/scBONITA2/scBONITA_output/rules_output/${DATASET_NAME}_${NUM_GENES}_rules/${DATASET_NAME}_${NUM_GENES}.graphml_${DATASET_NAME}_${NUM_GENES}_rules.txt" \
        --simulated_rules "/home/emoeller/github/scBONITA2/network_simulation/data/${DATASET_NAME}_${NUM_GENES}_rules.txt" \
        --data_file "/home/emoeller/github/scBONITA2/input/${DATASET_NAME}_${NUM_GENES}.csv" \
        --output_file "${DATASET_NAME}_${NUM_GENES}.log"

    echo ""
done
