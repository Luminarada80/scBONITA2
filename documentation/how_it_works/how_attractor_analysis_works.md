# Attractor Analysis

The attractor analysis portion of scBONITA2 consists of two main steps: **simulating individual cells** and **clustering the cells** based on the similarity of their signaling trajectories for the network.

## Simulating cell trajectories

First, the network pickle file containing the information about the pathway model is loaded. The network object stores the data for each gene in the network, the information about each node in the network, and the name of the network.

The cell trajectories are simulated by randomly choosing `num_simulations` columns in the dataset. For each chosen column (cell), a vectorized_run_simulation function similar to the one used for calculating importance scores is used to simulate the Boolean network model, starting from the chosen cell's expression values for the genes in the network. This process is multithreaded to increase the computational efficiency when simulating a large number of cells. Once a cell trajectory has been simulated, the resulting trajectory is stored as a csv file.
