import networkx as nx

original_graph = nx.read_graphml('temp_scSEQ_antiviral_signaling_net.graphml')

# Create a new graph of the same type as the original (directed/undirected)
new_graph = type(original_graph)()

# Iterate over all nodes in the original graph
for node, data in original_graph.nodes(data=True):
    gene_symbol = data.get('gene_symbol')  # Get the gene_symbol for the node
    if gene_symbol:
        # Add a new node to the new graph using gene_symbol as the ID
        new_graph.add_node(gene_symbol, **data)  # Copy all attributes

# Iterate over all edges in the original graph and add them to the new graph
for u, v, data in original_graph.edges(data=True):
    gene_symbol_u = original_graph.nodes[u].get('gene_symbol')
    gene_symbol_v = original_graph.nodes[v].get('gene_symbol')
    if gene_symbol_u and gene_symbol_v:
        new_graph.add_edge(gene_symbol_u, gene_symbol_v, **data)

# At this point, new_graph contains nodes identified by their gene symbols
# Path to the new GraphML file
output_path = 'modified_network.graphml'

# Write the new graph to a GraphML file
nx.write_graphml(new_graph, output_path)
