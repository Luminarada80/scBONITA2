#!/bin/bash -l

DATASET_NAME="simulated_data"
NUM_GENES="25"

INPUT_DIR=""

python3 "/home/emoeller/github/scBONITA2/network_simulation/network_simulation.py" \
    --dataset_name "${DATASET_NAME}" \
    --num_genes "${NUM_GENES}"

# Run the pipeline class 
python3 "/home/emoeller/github/scBONITA2/code/pipeline_class.py" \
    --data_file "input/${DATASET_NAME}_${NUM_GENES}.csv" \
    --dataset_name "${DATASET_NAME}_${NUM_GENES}" \
    --datafile_sep "," \
    --network_files "${DATASET_NAME}_${NUM_GENES}.graphml" \
    --binarize_threshold 0.01 \
    --get_kegg_pathways "False" \
    --minimum_overlap 1


# Run the rule checker
python3 "/home/emoeller/github/scBONITA2/network_simulation/check_rules_output.py" \
    --scbonita_rules "/home/emoeller/github/scBONITA2/scBONITA_output/rules_output/${DATASET_NAME}_${NUM_GENES}_rules/${DATASET_NAME}_${NUM_GENES}.graphml_${DATASET_NAME}_${NUM_GENES}_rules.txt" \
    --simulated_rules "/home/emoeller/github/scBONITA2/network_simulation/data/${DATASET_NAME}_${NUM_GENES}_rules.txt" \
    --data_file "/home/emoeller/github/scBONITA2/input/${DATASET_NAME}_${NUM_GENES}.csv" \
    --output_file "${DATASET_NAME}_${NUM_GENES}.log" \
