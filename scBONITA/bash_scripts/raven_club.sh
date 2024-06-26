#!/bin/bash

# -------------- User Input --------------

# IMPORTANT!!! MAKE SURE THAT THERE ARE NO SPACES IN FILE NAMES

# Which parts do you want to run? Set True to run or False to skip
    # Rule determination must be run prior to importance score, importance score must be run prior to relative abundance
RUN_RULE_DETERMINATION=False
RUN_IMPORTANCE_SCORE=True
RUN_RELATIVE_ABUNDANCE=True
RUN_ATTRACTOR_ANALYSIS=True
RUN_CELL_MAPPING=True

# Set the directory path for where the main scBONITA files are found
HOME=/home/emoeller/github/scBONITA/

# General Arguments (Required for all steps)
DATA_FILE="$HOME/input/datasets/raven_data_covid/Club_Cells_data.csv"
DATASET_NAME="raven_club"
DATAFILE_SEP=","
#  "04010" "04370" "04630" "04668" "04066" "04020" "04151" "04150" "00010" "00020" "04060" "04512" "04514" "04670" "04625" "04062"  "04810"
KEGG_PATHWAYS=() # Enter KEGG pathway codes or leave blank to find all pathways with overlapping genes
FIND_PATHWAYS=False
CUSTOM_PATHWAYS=("05171_modified.graphml") #("modified_network.graphml") #Put custom networks in the scBONITA folder
BINARIZE_THRESHOLD=0.01 # Data points with values above this number will be set to 1, lower set to 0
ORGANISM_CODE="hsa" # Organism code in front of KEGG pathway numbers
#"05171_modified.graphml"

# Relative Abundance arguments
METADATA_FILE="$HOME/input/metadata_files/Club_Cells_meta.txt"
METADATA_SEP=" "
HEADER="y" # Does the metadata file contain a header before the entries start?
OVERWRITE="y" # Do you want to overwrite the files generated for each of your different experimental groups?
CELL_NAME_COL=1 # What column contains the cell names (first column = 0)
GROUP_INDICES=(2 3)

# Specify the control groups and experimental groups that you want to compare
    # 1st entry in control is compared to 1st entry in experimental, 2nd entry compared to 2nd entry, etc.
    # Separate multiple groups with "_" (e.g. if a cell has two groups CTRL and Adult, combine as CTRL_Adult)
        #CONTROL_GROUPS=("CTRL_Pediatric" "CTRL_Adult")
        #EXPERIMENTAL_GROUPS=("SARS-CoV-2_Pediatric" "SARS-CoV-2_Adult")

CONTROL_GROUPS=("CTRL_Pediatric" "CTRL_Adult")
EXPERIMENTAL_GROUPS=("SARS-CoV-2_Pediatric" "SARS-CoV-2_Adult")

# -------------- End of user input --------------

#  ----------------------------
# |     RULE DETERMINATION     |
#  ----------------------------

if [ "$RUN_RULE_DETERMINATION" = "True" ]; then
    echo "Running Rule Determination..."


    if [ ${#KEGG_PATHWAYS[@]} -gt 0 ]; then
        echo "Running with KEGG pathways"

        # Using a list of KEGG pathways:
        KEGG_PATHWAYS_ARGS="${KEGG_PATHWAYS[@]}"

        /home/emoeller/anaconda3/envs/scBonita/bin/python "$HOME"/scBONITA/pipeline_class.py \
            --data_file "$DATA_FILE" \
            --dataset_name "$DATASET_NAME" \
            --datafile_sep "$DATAFILE_SEP" \
            --list_of_kegg_pathways $KEGG_PATHWAYS_ARGS \
            --binarize_threshold $BINARIZE_THRESHOLD \
            --organism $ORGANISM_CODE

    elif [ "$FIND_PATHWAYS" = "True" ]; then
        echo "No KEGG pathways specified, finding kegg pathways with overlapping genes..."
        /home/emoeller/anaconda3/envs/scBonita/bin/python "$HOME"/scBONITA/pipeline_class.py \
        --data_file "$DATA_FILE" \
        --dataset_name "$DATASET_NAME" \
        --datafile_sep "$DATAFILE_SEP" \
        --get_kegg_pathways True \
        --binarize_threshold $BINARIZE_THRESHOLD \
        --organism $ORGANISM_CODE
    
    else
        echo "Find pathways False and no KEGG pathways specified, checking for custom pathways"
    fi

    # Using a custom network saved to the scBONITA directory:

    # Check and execute for Custom Pathways if the array is not empty
    if [ ${#CUSTOM_PATHWAYS[@]} -gt 0 ]; then
        echo "Running with Custom Pathways..."
        
        CUSTOM_PATHWAYS_ARGS=""
        for pathway in "${CUSTOM_PATHWAYS[@]}"; do
            CUSTOM_PATHWAYS_ARGS+="--network_files $pathway "
        done

        /home/emoeller/anaconda3/envs/scBonita/bin/python "$HOME"/scBONITA/pipeline_class.py \
        --data_file "$DATA_FILE" \
        --dataset_name "$DATASET_NAME" \
        --datafile_sep "$DATAFILE_SEP" \
        $CUSTOM_PATHWAYS_ARGS \
        --binarize_threshold $BINARIZE_THRESHOLD \
        --get_kegg_pathways "False"
    else
        echo "No Custom Pathways specified, skipping this part..."
    fi
fi

#  --------------------------------------
# |     IMPORTANCE SCORE CALCULATION     |
#  --------------------------------------

if [ "$RUN_IMPORTANCE_SCORE" = "True" ]; then
    echo "Running Importance Score Calculation..."

    if [ ${#KEGG_PATHWAYS[@]} -gt 0 ]; then
        echo "Running with KEGG pathways"
        KEGG_PATHWAYS_ARGS="${KEGG_PATHWAYS[@]}"
    else
        echo "No KEGG Pathways specified"
    fi

    if [ ${#CUSTOM_PATHWAYS[@]} -gt 0 ]; then
        echo "Running with Custom Pathways..."
        CUSTOM_PATHWAYS_ARGS=""
        for pathway in "${CUSTOM_PATHWAYS[@]}"; do
            CUSTOM_PATHWAYS_ARGS+="--list_of_kegg_pathways $pathway "
        done
    else
        echo "No Custom Pathways specified"
    fi

    # Build the command with common arguments
    CMD="/home/emoeller/anaconda3/envs/scBonita/bin/python \"$HOME/scBONITA/importance_scores.py\" \
        --dataset_name \"$DATASET_NAME\""

    # Add KEGG pathways if available
    if [ -n "$KEGG_PATHWAYS_ARGS" ]; then
        CMD+=" --list_of_kegg_pathways $KEGG_PATHWAYS_ARGS"
    fi

    # Add custom pathways if available
    if [ -n "$CUSTOM_PATHWAYS_ARGS" ]; then
        CMD+=" $CUSTOM_PATHWAYS_ARGS"
    fi

    # Execute the command
    eval $CMD
fi

#  -----------------------------------------
# |     RELATIVE ABUNDANCE CALCULATIONS     |
#  -----------------------------------------

if [ "$RUN_RELATIVE_ABUNDANCE" = "True" ]; then
    echo "Running Relative Abundance Calculations..."

    GROUP_INDICES_ARGS="${GROUP_INDICES[@]}"

    # Check that both arrays have the same length
    if [ ${#CONTROL_GROUPS[@]} -ne ${#EXPERIMENTAL_GROUPS[@]} ]; then
        echo "Control and Experimental groups arrays do not match in length!"
        exit 1
    fi

    if [ ${#KEGG_PATHWAYS[@]} -gt 0 ]; then
        echo "Running with KEGG pathways"
        KEGG_PATHWAYS_ARGS="${KEGG_PATHWAYS[@]}"
    else
        echo "No KEGG Pathways specified"
    fi

    if [ ${#CUSTOM_PATHWAYS[@]} -gt 0 ]; then
        echo "Running with Custom Pathways..."
        CUSTOM_PATHWAYS_ARGS=""
        for pathway in "${CUSTOM_PATHWAYS[@]}"; do
            CUSTOM_PATHWAYS_ARGS+="--list_of_kegg_pathways $pathway "
        done
    fi

    # Loop through the control and experimental groups
    for (( i=0; i<${#CONTROL_GROUPS[@]}; i++ )); do

        # Extract the current pair of control and experimental group
        CONTROL_GROUP=${CONTROL_GROUPS[$i]}
        EXPERIMENTAL_GROUP=${EXPERIMENTAL_GROUPS[$i]}

        # Build the command with common arguments
        CMD="/home/emoeller/anaconda3/envs/scBonita/bin/python \"$HOME/scBONITA/relative_abundance.py\" \
            --dataset_name \"$DATASET_NAME\" \
            --dataset_file \"$DATA_FILE\" \
            --metadata_file \"$METADATA_FILE\" \
            --metadata_sep \"$METADATA_SEP\" \
            --dataset_sep \"$DATAFILE_SEP\" \
            --control_group \"$CONTROL_GROUP\" \
            --experimental_group \"$EXPERIMENTAL_GROUP\" \
            --cell_name_index $CELL_NAME_COL \
            --group_indices $GROUP_INDICES_ARGS \
            --header \"$HEADER\" \
            --overwrite \"$OVERWRITE\" \
            --organism \"$ORGANISM_CODE\""

        # Add custom pathways if available
        if [ -n "$CUSTOM_PATHWAYS_ARGS" ]; then
            CMD+=" $CUSTOM_PATHWAYS_ARGS"
        elif [ -n "$KEGG_PATHWAYS_ARGS" ]; then
            CMD+=" --list_of_kegg_pathways $KEGG_PATHWAYS_ARGS"
        fi

        # Execute the command
        eval $CMD
    done
fi


#  --------------------------------------
# |          ATTRACTOR ANALYSIS          |
#  --------------------------------------

# Runs the attractor analysis, requires importance score calculations
if [ "$RUN_ATTRACTOR_ANALYSIS" = "True" ]; then
    echo "Running Attractor Analysis..."

    /home/emoeller/anaconda3/envs/scBonita/bin/python "$HOME"/scBONITA/attractor_analysis.py \
        --dataset_name "$DATASET_NAME"
fi

#  --------------------------------------
# |             CELL MAPPING             |
#  --------------------------------------

# Maps each cell to the attractors from each network that best matches its gene expression
# Requires attractor analysis
if [ "$RUN_CELL_MAPPING" = "True" ]; then
    echo "Running Cell Mapping..."

    /home/emoeller/anaconda3/envs/scBonita/bin/python "$HOME"/scBONITA/map_cells_to_attractor_clusters.py \
        --dataset_name "$DATASET_NAME"
fi
