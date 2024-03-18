22# scBONITA2
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

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/f44738a4-93d0-4266-8bd7-2f3376ffabbc)

3) The cells are clustered by randomly shuffling the columns of the dataset and binning the expression data for each gene based on whether the majority of the cells within the bin are 1 or 0.

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/15878e85-dbdf-4cd2-9eb3-aad8890ae2fb)

4) The genetic algorithm is ran to find a set of rules that has a low error. Each of the nodes in the network has other nodes signaling to it. These incoming nodes can be connected by Boolean rules, for example if node 1 and node 2 signal to node 3, the possibilities are that node 1 AND node 2 activate node 3 or that node 1 OR node 2 activates node 3.

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/d2970c94-294b-4906-b457-11b1190fe302)
![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/cd262972-b341-4d3f-966d-2a5f97e39fec)

5) Each of the nodes could also be inhibitory, for example if node 1 is inhibitory the rule might be NOT node 1 AND node 2 or NOT node 1 OR node 2. The possibilities for each node are generated and one rule for each node is randomly sampled to create individuals for the genetic algorithm, and the error for each one of the rules is calculated by seeing if the predicted logic fits the data. For example, if the predicted rule is node 1 AND node 2 activate node 3, and the data for cell X is node 1 = 0, node 2 = 1, node 3 = 1, then that would be an error because for an AND rule, we would expect that 0 AND 1 = 0.

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/1ab0cf31-f5a1-4d3f-9ec6-7c2b0a0326f2)

6) We repeat this for each cell cluster in the dataset, adding up the total error as the fitness of the individual. The higher the error, the lower the fitness. The population for the genetic algorithm is made of many individual random rulesets, and over time the individuals with the lowest error are selected.

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/8bab7b15-3c7a-41e4-8a06-d7a772a7325f)

9) Once the best individuals from the genetic algorithm are selected, a rule refinement method is performed on for each node in the best individuals, where the other rule possibilities are considered for nodes with a high error. This is done to optimize rules with high error that were not optimized during the genetic algorithm. The rulesets are formatted and output to the rules_output directory for the project.

### Importance Score Calculations
1) Once the ruleset is inferred, the flow of a signal through the network is simulated. We sample the cells in the dataset and set the state of the nodes in the network to the expression of the cell (if the cell state for a gene is 1, we set the initial state of the node to 1). We then apply the rules synchronously across the network and keep track of the state of each node in the network. The network is simulated until a cycle of states called an attractor is found. 

![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/e4863904-0f3c-406f-b22e-69acaddac292)

2) Each of the nodes is iteratively knocked out (expression is forced to be 0 throughout the simulation) or knocked in (expression is forced to be 1 throughout the simulation). The network is simulated again for the same set of cell states and the attractors recorded. Once the knock out and knock in simulations are conducted for each node in the network, the importance of each node is scored by recording the number of nodes that are altered in the attractor cycle when the nodes are knocked out and knocked in compared to normal. The greater the number of differences, the more important the node is. If a node changes the signaling pattern of the network greatly when its expression is altered, then under- or over-expression of that gene will have a greater impact on the pathway compared to a node that does not change the signaling pattern greatly.

   ![image](https://github.com/Luminarada80/scBONITA2/assets/140645994/1fb129d5-fd18-4293-beee-7afc571bfe2b)

> In this example, knocking out node 3 does not alter the expression as much as knocking in node 3


