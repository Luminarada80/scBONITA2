# scBONITA2
Infers Boolean molecular signaling networks using scRNAseq data and prior knowledge networks, performs attractor analysis, and calculates the importance of each node in the network

## Setup:

**Cloning repository**
Go the the directory where you want to install this project and enter `git clone https://github.com/Luminarada80/scBONITA2.git` and enter your username and password.

**Setting up your conda environment**
1) Install Anaconda from https://www.anaconda.com/download
2) Run the Anaconda installer in the terminal on Linux or WSL (Windows Subsystem for Linux)
   - `bash Anaconda3-20204.09-0-Linux-x86_64.sh` (if the file you downloaded is different, use that file name)
   - Follow the prompts to install
3) Once Anaconda is installed, close and re-open your terminal. You should see `(base)` before your username
4) To create the correct conda environment, navigate to the `scBONITA2` directory that you cloned from before in the terminal. Once at the correct directory, enter `conda env create --file spec-file.txt --name scBonita`
5) Once conda has finished working, you can confirm that the environment was created by entering `conda activate scBonita`. This will switch you into the correct environment to work with scBONITA.
   - A conda environment is basically a pre-packaged python installation that has a specific python version and package list that works with the code. This makes it so that you don't have to install each required package one-by-one, and you can have different package versions by having different conda environments

**Testing that scBONITA is working**
1) Ensure that the `scBonita` conda environment is active (or enter `conda activate scBonita` to activate)
2) Navigate to `scBONITA2/scBONITA` in your terminal
3) Run the test data for the scBONITA pipeline by entering `bash bash_scripts/local_george_hiv.sh`
   - The `local_george_hiv.sh` file can be copied and modified to run different datasets / conditions
   - When running scBONITA, place bash scripts into this folder to run. scBONITA will automatically create the necessary output files in a directory called `scBONITA_output`.
   - Make sure that you have the `george_data` directory downloaded 

## Running scBONITA
**BASH Script** 

scBONITA runs using a BASH script that allows you to specify all of your parameters in one place. These BASH scripts are found in the `scBONITA2/scBONITA/bash_scripts` directory. To run using your own data, simply copy the bash file and modify the values.  

Each of the steps for running scBONITA is modular, meaning that you can run the rule determination step first. This means that you dont have to run it all at the same time (although you do need to run it in order).

You have to run Rule Determination $\rightarrow$ Importance Score $\rightarrow$ Relative Abundance $\rightarrow$ Attractor Analysis, but not at the same time. I recommend running Rule Determination first to make sure your environment and paths are correct.


```bash
# IMPORTANT!!! MAKE SURE THAT THERE ARE NO SPACES IN FILE NAMES

# Which parts do you want to run? Set True to run or False to skip
    # Rule determination must be run prior to importance score, importance score must be run prior to relative abundance
RUN_RULE_DETERMINATION=True
RUN_IMPORTANCE_SCORE=False
RUN_RELATIVE_ABUNDANCE=False
RUN_ATTRACTOR_ANALYSIS=False

# General Arguments (Required for all steps)
DATA_FILE="../../george_data/hiv_dataset/HIV_dataset_normalized_integrated_counts.csv"
DATASET_NAME="george_hiv"
DATAFILE_SEP=","
KEGG_PATHWAYS=("04370") # Enter KEGG pathway codes or leave blank to find all pathways with overlapping genes. Separate like: ("hsa04670" "hsa05171")
CUSTOM_PATHWAYS=() #("modified_network.graphml") #Put custom networks in the input folder
BINARIZE_THRESHOLD=0.01 # Data points with values above this number will be set to 1, lower set to 0
MINIMUM_OVERLAP=1 # Specifies how many genes you want to ensure overlap with the genes in the KEGG pathways. Default is 25
ORGANISM_CODE="hsa" # Organism code in front of KEGG pathway numbers

# Relative Abundance arguments
METADATA_FILE="../../george_data/hiv_dataset/hiv_meta.txt"
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

For the relative abundance section, scBONITA expects a metadata file with the following format:

```
"0" "C1_Healthy_AAACGGGAGTCGTTTG.1" "Healthy"
"1" "C1_Healthy_AAACGGGTCAGGCGAA.1" "Healthy"
"2" "C1_Healthy_AAAGCAATCTCAACTT.1" "Healthy"
"3" "C1_Healthy_AACACGTAGAGCAATT.1" "Healthy"
"4" "C1_Healthy_AACCGCGCAATCCGAT.1" "Healthy"
"5" "C1_Healthy_AACCGCGTCCTATTCA.1" "Healthy"
```

You can specify the separator, cell name column, which columns have the name of the groups (the first column is 0), and whether or not there is a header on the first line. Specify the control groups and the experimental groups to compare for the relative abundance calculations.

## Understanding the Process:
### Network Processing:
1) scBONITA starts by finding either all KEGG pathways (if the `KEGG_PATHWAYS` variable in the bash file is blank) or the KEGG pathways specified in the `KEGG_PATHWAYS` variable in the bash file. It will download all pathway xml files into the `pathway_xml_files` directory specific to the organism to make parsing the pathways faster than streaming the data directly from the KEGG API.
2) Next, scBONITA will check that the dataset you provided has enough overlapping genes with the pathways specified by the user, or all KEGG pathways (default is >25 shared genes, you can modify this value in the BASH file if you need to).
3) Once the pathways are identified, the pathways are processed to create "_processed.graphml" files.

The output should look similar to this:

```
-----PARSING NETWORKS-----
        KEGG pathways = ['04370']
        Finding and formatting KEGG Pathways...
                Finding KEGG pathways...
                Parsing KEGG dict...
                        Reading in KEGG dictionary file...
                        Loaded KEGG code dictionary
                        Reading hsa dictionary file...
                Downloading any missing pathway xml files, this may take a while...
|████████████████████████████████████████| 359/359 [100%] in 0.2s (1892.88/s)
                Finding pathways with at least 1 genes that overlap with the dataset
|               Reading KEGG xml file    | ▁▃▅ 0/1 [0%] in 0s (0.0/s, eta: ?)
                        Pathway (0/1): 04370 Overlap: 30 Edges: 65
                Reading KEGG xml file
                        Pathway (0/1): 04370 Overlap: 59 Edges: 166
|████████████████████████████████████████| 1/1 [100%] in 0.1s (9.92/s)
                Adding graphml pathways to rule_inference object...
                Pathway: 04370 Overlap: 59 Edges: 166
                        Edges after processing: 158 Overlap: 59
```

### Loading the Data:
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

### Rule Inference
1) The network is loaded as a NetworkX object from the processed graphml file.
2) The cells are clustered by randomly shuffling the columns of the dataset and binning the expression data for each gene based on whether the majority of the cells within the bin are 1 or 0.

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/15878e85-dbdf-4cd2-9eb3-aad8890ae2fb)

3) The genetic algorithm runs to find a set of rules that has a low error. Each of the nodes in the network has other nodes signaling to it. These incoming nodes can be connected by Boolean rules, for example if node 1 and node 2 signal to node 3, the possibilities are that node 1 AND node 2 activate node 3 or that node 1 OR node 2 activates node 3.

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/d2970c94-294b-4906-b457-11b1190fe302)

4) Each of the nodes could also be inhibitory, for example if node 1 is inhibitory the rule might be NOT node 1 AND node 2 or NOT node 1 OR node 2. The possibilities for each node are generated and one rule for each node is randomly sampled to create individuals for the genetic algorithm, and the error for each one of the rules is calculated by seeing if the predicted logic fits the data. For example, if the predicted rule is node 1 AND node 2 activate node 3, and the data for cell X is node 1 = 0, node 2 = 1, node 3 = 1, then that would be an error because for an AND rule, we would expect that 0 AND 1 = 0.

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/1ab0cf31-f5a1-4d3f-9ec6-7c2b0a0326f2)

5) We repeat this for each cell cluster in the dataset, adding up the total error as the fitness of the individual. The higher the error, the lower the fitness. The population for the genetic algorithm is made of many individual random rulesets, and over time the individuals with the lowest error are selected.

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/8bab7b15-3c7a-41e4-8a06-d7a772a7325f)

6) Once the best individuals from the genetic algorithm are selected, a rule refinement method is performed on for each node in the best individuals, where the other rule possibilities are considered for nodes with a high error. This is done to optimize rules with high error that were not optimized during the genetic algorithm. The rulesets are formatted and output to the rules_output directory for the project.

### Importance Score Calculations
1) Once the ruleset is inferred, the flow of a signal through the network is simulated. We sample the cells in the dataset and set the state of the nodes in the network to the expression of the cell (if the cell state for a gene is 1, we set the initial state of the node to 1). We then apply the rules synchronously across the network and keep track of the state of each node in the network. The network is simulated until a cycle of states called an attractor is found. 

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/e4863904-0f3c-406f-b22e-69acaddac292)

2) Each of the nodes is iteratively knocked out (expression is forced to be 0 throughout the simulation) or knocked in (expression is forced to be 1 throughout the simulation). The network is simulated again for the same set of cell states and the attractors recorded. Once the knock out and knock in simulations are conducted for each node in the network, the importance of each node is scored by recording the number of nodes that are altered in the attractor cycle when the nodes are knocked out and knocked in compared to normal. The greater the number of differences, the more important the node is. If a node changes the signaling pattern of the network greatly when its expression is altered, then under- or over-expression of that gene will have a greater impact on the pathway compared to a node that does not change the signaling pattern greatly.

   ![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/1fb129d5-fd18-4293-beee-7afc571bfe2b)

> In this example, knocking out node 3 does not alter the expression as much as knocking in node 3


