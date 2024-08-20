# Rule Inference
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