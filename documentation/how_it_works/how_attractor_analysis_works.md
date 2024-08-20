# Attractor Analysis

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
