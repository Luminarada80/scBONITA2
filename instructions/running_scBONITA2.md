
# Running scBONITA:

1. scBONITA runs using a [bash script](bash_script_instructions.md) as a wrapper to pass user arguments into the Python script.

2. The [networks are downloaded and processed](how_it_works/how_network_processing_works.md) to create graphml files containing the genes that match your dataset, their connections to other genes, and whether that connection is activating or inhibiting.





## Loading the Data:
1) scBONITA first extracts the data, randomly sampling the cells if there are more than 15000 samples to limit the size of the dataset.
2) The matrix is then converted to a sparse matrix to condense the information, and the information is binarized to 1 or 0 based on the value of the `BINARIZE_THRESHOLD` variable in the bash file. Expression values above this threshold are set to 1 and values below this threshold are set to 0, indicating that the gene is either ON or OFF
3) The genes in the dataset are filtered to include variable genes with a coefficient of variation above the cv_threshold (default 0.001, this can be set as a pipeline.py argument in the bash file).
4) Once the dataset is loaded and configured, the processed network and data are used for determining the Boolean rules linking genes together.

Here is what a successful run looks like:

```
-----RULE INFERENCE-----
Pathway: 04370
Num nodes: 59

-----EXTRACTING AND FORMATTING DATA-----
Extracting cell expression data from "../../george_data/hiv_dataset/HIV_dataset_normalized_integrated_counts.csv"
        Loading all 3621 cells...
        Converting filtered data to numpy array...
        First 2 genes: ['PIK3CD', 'CASP9']
        First 2 cells: ['C1_Healthy_AAACGGGAGTCGTTTG.1', 'C1_Healthy_AAACGGGTCAGGCGAA.1']
        Number of genes: 59
        Number of cells: 3621
        Created sparse matrix
        Binarized sparse matrix
        Setting ruleset parameters
        Running rule inference for 04370
                Loading: hsa04370_processed.graphml
/home/emoeller/anaconda3/envs/scBonita/lib/python3.6/site-packages/scipy/stats/stats.py:4196: SpearmanRConstantInputWarning: An input array is constant; the correlation coefficent is not defined.
  warnings.warn(SpearmanRConstantInputWarning())
|████████████████████████████████████████| 59/59 [100%] in 0.0s (1838.83/s)
```

## Rule Inference
1) The network is loaded as a NetworkX object from the processed graphml file.
2) The cells are clustered by randomly shuffling the columns of the dataset and binning the expression data for each gene based on whether the majority of the cells within the bin are 1 or 0.
3) The genetic algorithm runs to find a set of rules that has a low error. Each of the nodes in the network has other nodes signaling to it. These incoming nodes can be connected by Boolean rules, for example if node 1 and node 2 signal to node 3, the possibilities are that node 1 AND node 2 activate node 3 or that node 1 OR node 2 activates node 3.

5) We repeat this for each cell cluster in the dataset, adding up the total error as the fitness of the individual. The higher the error, the lower the fitness. The population for the genetic algorithm is made of many individual random rulesets, and over time the individuals with the lowest error are selected.

6) Once the best individuals from the genetic algorithm are selected, a rule refinement method is performed on for each node in the best individuals, where the other rule possibilities are considered for nodes with a high error. This is done to optimize rules with high error that were not optimized during the genetic algorithm. The rulesets are formatted and output to the rules_output directory for the project.

Here is what the output from a successful rule output should look like:

```
-----CHUNKING DATASET-----
        Original Data Shape: 59 rows, 3621 columns
        Chunked Data Shape: 59 rows, 100 columns
        Coarse Chunked Data Shape: 59 rows, 50 columns

-----GENETIC ALGORITHM-----
ngen    nevals  avg     std     min     max
0       20      0.233   0.011   0.213   0.256
1       20      0.234   0.009   0.218   0.252
2       20      0.232   0.012   0.218   0.252
3       20      0.231   0.011   0.221   0.252
4       20      0.233   0.013   0.211   0.252

-----RULE REFINEMENT-----
Equivalent Ruleset 1 / 1
|████████████████████████████████████████| 59/59 [100%] in 2.6s (22.41/s)

Equivalent Ruleset 1
CDC42 = KDR
KDR = VEGFA
PTGS2 = not NFATC2
HSPB1 = MAPKAPK2 and MAPKAPK3
NFATC2 = (PPP3CA or PPP3CC) and PPP3CB
PTK2 = KDR
PXN = KDR
SH2D2A = (PLCG1 or KDR) and PLCG2
SRC = KDR
SHC2 = KDR
VEGFA = VEGFA
RAF1 = NRAS or HRAS and PRKCA
NOS3 = (PLCG2 or AKT1) and PLCG1
CASP9 = (AKT1 or AKT3) and AKT2
.
.
.
MAP2K1 = RAF1
MAP2K2 = RAF1
HRAS = SPHK2 and SPHK1
KRAS = SPHK1 and SPHK2
NRAS = SPHK2 and SPHK1
RAC1 = PIK3R2 or PIK3CB and PIK3R1
RAC2 = (PIK3R2 or PIK3CD) and PIK3R1
RAC3 = (PIK3R2 or PIK3CD) and PIK3R1
Refined Error:

        Average = 0.08474576271186439
        Stdev = 0.14517950758075673
        Max = 0.57
        Min = 0.0
        Percent_self_loops = 2.0%

Rule inference complete, saving to george_hiv_hsa04370.ruleset.pickle
``` 

## Importance Score Calculations
1) Once the ruleset is inferred, the flow of a signal through the network is simulated. We sample the cells in the dataset and set the state of the nodes in the network to the expression of the cell (if the cell state for a gene is 1, we set the initial state of the node to 1). We then apply the rules synchronously across the network and keep track of the state of each node in the network. The network is simulated until a cycle of states called an attractor is found. 

2) Each of the nodes is iteratively knocked out (expression is forced to be 0 throughout the simulation) or knocked in (expression is forced to be 1 throughout the simulation). The network is simulated again for the same set of cell states and the attractors recorded. 

3) Once the knock out and knock in simulations are conducted for each node in the network, the importance of each node is scored by recording the number of nodes that are altered in the attractor cycle when the nodes are knocked out and knocked in compared to normal. 
    - The greater the number of differences, the more important the node is. 
    - If a node changes the signaling pattern of the network greatly when its expression is altered, then under- or over-expression of that gene will have a greater impact on the pathway compared to a node that does not change the signaling pattern greatly.

Here is what the output from a successful importance score calculation step should look like:

```
 --------------------------------------------------
|     RUNNING IMPORTANCE SCORES FOR GEORGE_HIV     |
 --------------------------------------------------

Loading: george_hiv_hsa04370.ruleset.pickle
Calculating importance score for network hsa04370

-----RUNNING NETWORK SIMULATION-----
|████████████████████████████████████████| 59/59 [100%] in 22.8s (2.59/s)

-----CALCULATING IMPORTANCE SCORES-----
|████████████████████████████████████████| 59/59 [100%] in 0.9s (68.14/s)
Saving importance scores to file: 04370_george_hiv_importance_scores.txt
Saving importance score figures
Saving network object as a pickle file
```

## Relative Abundance

**User input in the BASH file**
```bash
# Relative Abundance arguments
METADATA_FILE="../input/george_data/hiv_dataset/hiv_meta.txt"
METADATA_SEP=" "
HEADER="n" # Does the metadata file contain a header before the entries start?
OVERWRITE="n" # Do you want to overwrite the files generated for each of your different experimental groups?
CELL_NAME_COL=1 # What column contains the cell names (first column = 0)
GROUP_INDICES=(2)

# Specify the control groups and experimental groups that you want to compare
    # 1st entry in control is compared to 1st entry in experimental, 2nd entry compared to 2nd entry, etc.
CONTROL_GROUPS=("Healthy")
EXPERIMENTAL_GROUPS=("HIV")
```

**Metadata File**

```
"0" "C1_Healthy_AAACGGGAGTCGTTTG.1" "Healthy"
"1" "C1_Healthy_AAACGGGTCAGGCGAA.1" "Healthy"
"2" "C1_Healthy_AAAGCAATCTCAACTT.1" "Healthy"
"3" "C1_Healthy_AACACGTAGAGCAATT.1" "Healthy"
"4" "C1_Healthy_AACCGCGCAATCCGAT.1" "Healthy"
```


**Argument Explanation:** 

Here, the dataset is being split based on the metadata file `hiv_meta.txt`. Columns in the file are separated with a space (`" "`), there is no header, and the names of the cells are in column 1 (indexing starts at 0). The group each cell belongs to is found in column 2, and the two groups of interest are "Healthy" and "HIV".

1) First, the columns in the dataset are split into the two groups "Healthy" and "HIV" and saved to separate csv files in the same location as the original data file. If you have already created the split group data files in the past and want to overwrite them, specify `OVERWRITE="y"` in the BASH file.

```
 ----------------------------------------------------------
|     RELATIVE ABUNDANCE FOR GEORGE_HIV HEALTHY VS HIV     |
 ----------------------------------------------------------

----- Splitting Data Files -----
        Group 1: Healthy
        Group 2: HIV

----- Saving Group Data Files -----
                Writing data to the group file...
                Writing data to the group file...
                Using existing group network hsa04370 file for Healthy
                Using existing group network hsa04370 file for HIV
                Control Group: Healthy
                Experimental Group: HIV
```

2) Next, the networks of interest are loaded for both groups

```
----- Loading george_hiv networks -----
        Loading CONTROL group networks
                Loaded network hsa04370_Healthy

        Loading EXPERIMENTAL group networks
                Loaded network hsa04370_HIV
```

3) Then, the relative abundance of each group is calculated along with a pathway modulation score. The pathway modulation score weighs the importance score and relative abundance of each node to determine if the pathway is significantly altered. It uses bootstrapping to determine if the pathway modulation score is significantly different compared to if the importance scores and relative abundances were randomly distributed among the genes.

> The distribution generated by bootstrapping can be viewed in the relative abundnace output folder

```
----- Calculating Relative Abundance between groups HIV and Healthy -----
        Network: hsa04370
                Pathway Modulation Score: 0.06087929272259514
                Calculating p-value with bootstrapping:
                        P-value: 0.5676
                        -log10(P-value): 0.2459576131380494
```

## Attractor Analysis

> NOTE: This is a current work in progress, I am switching to a different method to map cells to attractors.

This process attempts to find attractors that are similar to one another and create a representative attractor that best represents the different signaling states of each cell. Each cell is compared to each attractor to find the attractor that best represents that cell. Attractors that do not match most closely to at least one cell are removed. Hamming distance is then used to compare each attractor to generate a distance matrix. The attractors are clustered using Hierarchical clustering. Each cell is then assigned to a cluster based on which cluster contains the attractor it is most similar to. The attractor that best matches the greatest number of cells is used as the representative attractor for that cluster.

Here is what a successful run should look like:

```
----- ATTRACTOR ANALYSIS -----

Network: hsa04370
        Generating attractors...
        Calculating hamming distance between cells and attractors
        Transposed dataset shape: (59, 500)
                Extra dimension in dataset, squeezing...
                Transposed dataset shape: (500, 59)
                Nodes: 59
                Cells: 500
                Attractors: 500
        Generating attractor distance matrix...
        Clustering the attractors...
                Clustering cutoff value = 14.75 (<=20% * number of genes 59)
        Calculating Hamming distance between cells and clustered attractors
Number of cells in the full dataset: 3621

-----PATHWAY ANALYSIS RESULTS -----

NETWORK HSA04370
        Attractor 1 contains 3619 cells (99.945%)
        Attractor 2 contains 2 cells (0.055%)
        Saved representative attractor 1
        Saved representative attractor 2

        Saved attractor analysis results to "attractor_analysis_output/george_hiv_attractors/hsa04370_attractors

Adding representative attractor map to network pickle files:
        File: george_hiv_hsa04370.network.pickle
```

## Dynamic Time Warping

First, run the `simulate_cell.py` script in the `scBONITA2/scBONITA` directory to generate cell trajectories for the networks of interest

`python3 simulate_cell.py --dataset_name george_hiv --network_name hsa04370 --num_cells 50`

Once the cell trajectories are created, we can use the `dynamic_time_warping.py` script to group 
the cells by how similar their trajectories are.

`python3 dynamic_time_warping.py --dataset_name george_hiv --network_name hsa04370`

This code first uses dynamic time warping to compute how different each cell's trajectory is from
every other cell. Once the pairwise differences are calculated, the cells are clustered using
hierarchical clustering to group the cells. The user will be prompted to select the number of 
clusters to split the cells into. Finally, a summary of the average signaling for each gene at each
simulation step is generated to show the average trajectory for that cluster.