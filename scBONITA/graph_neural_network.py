import networkx as nx
import pandas as pd
import torch
from torch_geometric.data import Data
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from torch_geometric.nn import GATConv
from torch_geometric.loader import DataLoader
from argparse import ArgumentParser
import pickle
import os
import random
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean
from file_paths import file_paths
import statistics
from alive_progress import alive_bar

class GNN(torch.nn.Module):
    def __init__(self, num_node_features, hidden_channels, output_dim):
        super(GNN, self).__init__()
        self.conv1 = GCNConv(num_node_features, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        self.conv3 = GCNConv(hidden_channels, hidden_channels)
        self.dropout = torch.nn.Dropout(p=0.5)
        self.linear = torch.nn.Linear(hidden_channels, output_dim)
        self.residual = torch.nn.Linear(num_node_features, hidden_channels)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x_residual = self.residual(x)
        x = self.conv1(x, edge_index)
        x = F.relu(x + x_residual)
        x = self.dropout(x)
        x = self.conv2(x, edge_index)
        x = F.relu(x)
        x = self.dropout(x)
        x = self.conv3(x, edge_index)
        x = F.relu(x)
        x = self.dropout(x)
        x = self.linear(x)
        return x

class GAT(torch.nn.Module):
    def __init__(self, num_node_features, hidden_channels, output_dim):
        super(GAT, self).__init__()
        self.conv1 = GATConv(num_node_features, hidden_channels)
        self.conv2 = GATConv(hidden_channels, hidden_channels)
        self.conv3 = GATConv(hidden_channels, hidden_channels)
        self.dropout = torch.nn.Dropout(p=0.5)
        self.linear = torch.nn.Linear(hidden_channels, output_dim)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.dropout(x)
        x = self.conv2(x, edge_index)
        x = F.relu(x)
        x = self.dropout(x)
        x = self.conv3(x, edge_index)
        x = F.relu(x)
        x = self.dropout(x)
        x = self.linear(x)
        return x

class EarlyStopping:
    def __init__(self, patience=10, min_delta=0):
        self.patience = patience
        self.min_delta = min_delta
        self.best_loss = None
        self.counter = 0
        self.early_stop = False

    def __call__(self, val_loss):
        if self.best_loss is None:
            self.best_loss = val_loss
        elif val_loss < self.best_loss - self.min_delta:
            self.best_loss = val_loss
            self.counter = 0
        else:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True

def calculate_dtw(actual, predicted):
    distance, path = fastdtw(actual, predicted, dist=euclidean)
    return distance
    

def prepare_data(network, starting_states, trajectories):
    print(f'Preparing data for the GNN')
    edge_index = torch.tensor(list(network.edges)).t().contiguous()
    data_list = []

    for key in starting_states:
        start_state = starting_states[key].values.astype(float)
        trajectory = trajectories[key].values.astype(float)
        
        x = torch.tensor(start_state, dtype=torch.float).unsqueeze(1)
        y = torch.tensor(trajectory, dtype=torch.float)
        
        data = Data(x=x, edge_index=edge_index, y=y)
        data_list.append(data)

    return data_list

def train_model(data_list, model, epochs=100, lr=0.01, weight_decay=1e-5, patience=10):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=20, gamma=0.1)
    loss_fn = torch.nn.BCEWithLogitsLoss()
    early_stopping = EarlyStopping(patience=patience)

    # Lists to store losses
    train_losses = []
    val_losses = []

    # Split data into training and validation sets
    train_data, val_data = train_test_split(data_list, test_size=0.2, random_state=42)

    for epoch in range(epochs):
        total_loss = 0
        model.train()  # Set model to training mode
        for data in train_data:
            optimizer.zero_grad()
            out = model(data)
            loss = loss_fn(out, data.y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        avg_train_loss = total_loss / len(train_data)
        train_losses.append(avg_train_loss)

        # Validation
        model.eval()  # Set model to evaluation mode
        total_val_loss = 0
        with torch.no_grad():
            for data in val_data:
                out = model(data)
                loss = loss_fn(out, data.y)
                total_val_loss += loss.item()
        avg_val_loss = total_val_loss / len(val_data)
        val_losses.append(avg_val_loss)

        # Step the learning rate scheduler
        scheduler.step()

        print(f'Epoch {epoch+1}, Train Loss: {avg_train_loss:.6f}, Validation Loss: {avg_val_loss:.6f}')

        early_stopping(avg_val_loss)
        if early_stopping.early_stop:
            print(f"Early stopping at epoch {epoch+1}")
            break

    return train_losses, val_losses

def read_data_from_directory(directory):
    print(f'Reading data from the files')
    starting_states = {}
    trajectories = {}
    gene_to_idx = {}
    cell_name_list = []
    idx = 0
    num_files = 0

    # First pass to build the gene_to_idx mapping
    for filename in os.listdir(directory):
        if filename.endswith("_trajectory.csv"):
            num_files += 1
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath, header=None)
            for gene in df.iloc[:, 0]:
                if gene not in gene_to_idx:
                    gene_to_idx[gene] = idx
                    idx += 1

    print(f'Found {num_files} trajectory files')
    
    # Second pass to convert gene names to indices and extract data
    for filename in os.listdir(directory):
        if filename.endswith("_trajectory.csv"):
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath, header=None)
            df.columns = ['Gene'] + [f'Time{i}' for i in range(df.shape[1] - 1)]
            df['Gene'] = df['Gene'].map(gene_to_idx)
            df.set_index('Gene', inplace=True)

            # Convert starting state and trajectory to integers
            df = df.apply(pd.to_numeric, errors='coerce')
            
            # Extract starting state and trajectory
            starting_state = df.iloc[:, 0]
            trajectory = df.iloc[:, 1:]
            
            cell_number = '_'.join(filename.split('_')[0:2])

            cell_name_list.append(cell_number)
            starting_states[filename] = starting_state
            trajectories[filename] = trajectory

            # if len(cell_name_list) < 5:
            #     print(f'{cell_number}')
            #     print(f'\t{starting_states[filename]}')
            #     print(f'\t{trajectories[filename]}')
    
    return cell_name_list, starting_states, trajectories, gene_to_idx

def min_max_scale(data):
    min_val = np.min(data)
    max_val = np.max(data)
    scaled_data = (data - min_val) / (max_val - min_val)
    return scaled_data

if __name__ == '__main__':
    parser = ArgumentParser()

    parser.add_argument(
        '--dataset_name',
        type=str,
        required=True,
        help='Number of genes to generate'
    )

    parser.add_argument(
        '--network_name',
        type=str,
        required=True,
        help='Number of cells to generate'
    )

    parser.add_argument(
        '--learn',
        type=str,
        required=True,
        help='Number of cells to generate'
    )

    results = parser.parse_args()

    dataset_name = getattr(results, 'dataset_name')
    network_name = getattr(results, 'network_name')
    learn = getattr(results, 'learn')

    # Directory containing the files
    directory = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/text_files'

    
    # Specifies the path to the correct network pickle file
    network_pickle_file = f'{file_paths["pickle_files"]}/{dataset_name}_pickle_files/network_pickle_files/{dataset_name}_{network_name}.network.pickle'

    # Read in the network pickle file from the path
    network_object = pickle.load(open(network_pickle_file, 'rb'))

    network = network_object.network

    edges = network.edges

    # Initialize the network
    new_network = nx.DiGraph()
    
    # Read in the trajectories as a tuple of starting state and resulting trajectory
    cell_name_list, starting_states, trajectories, gene_to_idx = read_data_from_directory(directory)

    # Convert gene names to indices for the edges
    edges = [(int(gene_to_idx[src]), int(gene_to_idx[dst])) for src, dst in edges]
    new_network.add_edges_from(edges)

    data_list = prepare_data(new_network, starting_states, trajectories)

    # Initialize model with output_dim matching the trajectory length
    output_dim = trajectories[next(iter(trajectories))].shape[1]  # Assuming all trajectories have the same length
    # model = GNN(num_node_features=1, hidden_channels=64, output_dim=output_dim)

    # Alternatively, use GAT model
    model = GAT(num_node_features=1, hidden_channels=64, output_dim=output_dim)

    if learn == "True":
        # Train the model with validation and early stopping, record losses
        train_losses, val_losses = train_model(data_list, model)

        # Save the trained model (optional, if you want to save and load later)
        torch.save(model.state_dict(), f'gat_{dataset_name}_{network_name}_model.pth')

        # Plot the training and validation loss
        plt.figure(figsize=(10, 6))
        plt.plot(train_losses, label='Training Loss')
        plt.plot(val_losses, label='Validation Loss')
        plt.xlabel('Epochs')
        plt.ylabel('Loss')
        plt.ylim((0,1))
        plt.title('Training and Validation Loss over Epochs')
        plt.legend()
        plt.show()

        plt.close()

    else:
        dtw_matrix = []

        model.load_state_dict(torch.load(f'gat_{dataset_name}_{network_name}_model.pth'))

        png_directory = f'{file_paths["trajectories"]}/{dataset_name}_{network_name}/png_files'

        # Set the model to evaluation mode
        model.eval()

        total_dtw_distances = []

        with alive_bar(len(data_list)) as bar:
            for i, sample_data in enumerate(data_list):
                # print(sample_data.x)
                # Make predictions
                with torch.no_grad():
                    output = model(sample_data)

                # Output is a tensor; convert it to a numpy array or list if needed
                output_array = output.numpy()

                # For binary classification, apply a sigmoid function to get probabilities
                predictions = torch.sigmoid(output).numpy()

                # Scale predictions between 0 and 1
                scaled_predictions = min_max_scale(predictions)

                # Calculate DTW between predictions and actual trajectory
                actual_trajectory = sample_data.y.numpy()
                
                # For each gene, determine if it is activated (1) or not (0) based on a threshold
                threshold = 0.5
                binary_predictions = (scaled_predictions > threshold).astype(int)

                dtw_distance = calculate_dtw(actual_trajectory, predictions)
                total_dtw_distances.append(dtw_distance)

                # Combine the starting state with the predicted trajectories
                starting_state = sample_data.x.numpy().reshape(-1, 1)  # Convert to 2D array if needed

                # Non-binarized combined data
                combined_data = np.concatenate((starting_state, scaled_predictions), axis=1)

                # Binarized combined data
                combined_data = np.concatenate((starting_state, binary_predictions), axis=1)

                # Create a heatmap
                plot = plt.figure(figsize=(12, 12))
                sns.heatmap(combined_data, cmap='Greys', yticklabels=gene_to_idx, xticklabels=True)
                plt.title(f'GNN Prediction for {cell_name_list[i]}')
                plt.xlabel('Time Steps')
                plt.ylabel('Genes')
                plt.xticks(fontsize=8)
                plt.yticks(fontsize=8)
                plt.tight_layout()
                # plt.show()

                # Saves a png of the results
                plot.savefig(f'{png_directory}/{cell_name_list[i]}_GAT_prediction.png', format='png')
                plt.close(plot)

                bar()

        average_dtw_distance = statistics.mean(total_dtw_distances)
        stdev_dtw_distance = statistics.stdev(total_dtw_distances)
        min_dtw_distance = min(total_dtw_distances)
        max_dtw_distance = max(total_dtw_distances)

        print(f'DTW distance btw prediction and trajectory: {average_dtw_distance}')
        print(f'Avg = {average_dtw_distance}')
        print(f'Stdev = {stdev_dtw_distance}')
        print(f'Min = {min_dtw_distance}')
        print(f'Max = {max_dtw_distance}')