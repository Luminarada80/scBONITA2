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
3) Run the test data for the scBONITA pipeline by entering `bash local_george_hiv.sh`
   - The `local_george_hiv.sh` file can be copied and modified to run different datasets / conditions

## Understanding the Process:
### Network Processing:
1) scBONITA starts by finding either all KEGG pathways (if the `KEGG_PATHWAYS` variable in the bash file is blank) or the KEGG pathways specified in the `KEGG_PATHWAYS` variable in the bash file. It will download all pathway xml files into the `pathway_xml_files` directory specific to the organism to make parsing the pathways faster than streaming the data directly from the KEGG API.
2) Next, scBONITA will check that the dataset you provided has enough overlapping genes with the pathways specified by the user, or all KEGG pathways (default is >25 shared genes).
3) Once the pathways are identified, the pathways are processed to create processed graphml files.

### Loading the Data:
1) scBONITA first extracts the data, randomly sampling the cells if there are more than 15000 samples to limit the size of the dataset.
2) The matrix is then converted to a sparse matrix to condense the information, and the information is binarized to 1 or 0 based on the value of the `BINARIZE_THRESHOLD` variable in the bash file. Expression values above this threshold are set to 1 and values below this threshold are set to 0, indicating that the gene is either ON or OFF
3) The genes in the dataset are filtered to include variable genes with a coefficient of variation above the cv_threshold (default 0.001, this can be set as a pipeline.py argument in the bash file).
4) Once the dataset is loaded and configured, the processed network and data are used for determining the Boolean rules linking genes together.

### Rule Inference
1) The network is loaded as a NetworkX object from the processed graphml file.
2) The cells are clustered by randomly shuffling the columns of the dataset and binning the expression data for each gene based on whether the majority of the cells within the bin are 1 or 0.
3) The genetic algorithm is ran to find a set of rules that has a low error. Each of the nodes in the network has other nodes signaling to it. These incoming nodes can be connected by Boolean rules, for example if node 1 and node 2 signal to node 3, the possibilities are that node 1 AND node 2 activate node 3 or that node 1 OR node 2 activates node 3. Each of the nodes could also be inhibitory, for example if node 1 is inhibitory the rule might be NOT node 1 AND node 2 or NOT node 1 OR node 2. The possibilities for each node are generated and one rule for each node is randomly sampled to create individuals for the genetic algorithm, and the error for each one of the rules is calculated by seeing if the predicted logic fits the data. For example, if the predicted rule is node 1 AND node 2 activate node 3, and the data for cell X is node 1 = 0, node 2 = 1, node 3 = 1, then that would be an error because for an AND rule, we would expect that 0 AND 1 = 0. We repeat this for each cell cluster in the dataset, adding up the total error as the fitness of the individual. The higher the error, the lower the fitness. The population for the genetic algorithm is made of many individual random rulesets, and over time the individuals with the lowest error are selected.
4) Once the best individuals from the genetic algorithm are selected, a rule refinement method is performed on for each node in the best individuals, where the other rule possibilities are considered for nodes with a high error. This is done to optimize rules with high error that were not optimized during the genetic algorithm. The rulesets are formatted and output to the rules_output directory for the project.

### Importance Score Calculations

