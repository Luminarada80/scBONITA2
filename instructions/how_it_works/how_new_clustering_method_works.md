# Dynamic Time Warping

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