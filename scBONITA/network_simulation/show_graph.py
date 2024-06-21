import networkx as nx
import matplotlib.pyplot as plt

def plot_graph_from_graphml(file_path):
    # Load the graph from a GraphML file
    G = nx.read_graphml(file_path)

    # Draw the graph
    plt.figure(figsize=(12, 8))  # Set the size of the plot
    nx.draw(G, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500, font_size=10)
    plt.title("Network View of the GraphML File")
    plt.show()

# Replace 'path_to_graphml_file.graphml' with the path to your GraphML file
plot_graph_from_graphml('boolean_network_conditions/random_boolean_network.graphml')
