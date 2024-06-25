#!/bin/bash

# Set the directory path for where the main scBONITA files are found
HOME=/home/emoeller/github/scBONITA

NUM_GENES=60
NUM_CELLS=1
NUM_CHUNKS=2000
ALLOW_MISMATCHES=false

until [ $NUM_GENES -gt 100 ]; do

    TOTAL_CELLS=$(($NUM_CELLS * $NUM_CHUNKS))

    echo "Running network simulation with ${NUM_GENES} genes and ${TOTAL_CELLS} cells"

    # Generate a network for the current number of genes
    /home/emoeller/anaconda3/envs/scBonita/bin/python "$HOME"/scBONITA/network_simulation/network_simulation.py \
    --num_genes $NUM_GENES \
    --num_cells $NUM_CELLS \
    --num_chunks $NUM_CHUNKS \
    --allow_mismatches "False"

    # Arguments for scBONITA
    DATA_FILE="$HOME/input/datasets/sim_dataset_${NUM_GENES}g_${TOTAL_CELLS}c.csv"
    DATASET_NAME="sim_data_${NUM_GENES}g_${TOTAL_CELLS}c"
    CUSTOM_PATHWAYS="$HOME/input/custom_graphml_files/sim_network_${NUM_GENES}g_${TOTAL_CELLS}c.graphml"
    BINARIZE_THRESHOLD=0.01

    # Run scBONITA rule inference on the dataset with the current number of genes
    /home/emoeller/anaconda3/envs/scBonita/bin/python "$HOME"/scBONITA/pipeline_class.py \
    --data_file "$DATA_FILE" \
    --dataset_name "$DATASET_NAME" \
    --datafile_sep "," \
    --network_files "$CUSTOM_PATHWAYS" \
    --binarize_threshold $BINARIZE_THRESHOLD \
    --get_kegg_pathways "False"
    
    # Increase the number of genes by 10
    NUM_GENES=$(($NUM_GENES + 10))
    
done

echo "Processing results using check_rules_output.py"

/home/emoeller/anaconda3/envs/scBonita/bin/python "$HOME"/scBONITA/network_simulation/check_rules_output.py