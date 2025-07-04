from pyrbd_plusplus import *

if __name__ == "__main__":
    # Example usage of minimalcuts_optimized
    directory = "topologies/Germany_17"  # Replace with your directory
    topo = "Germany_17"  # Replace with your topological name

    # Call the function to evaluate minimal cuts
    G, pos, lable = read_graph(directory, topo)
    
    pathsets = minimalpaths(G, 0, 1)
    
    print("Pathsets:", pathsets)
    
    # Read the graph
    directory = "topologies/Nobel_EU"  # Replace with your directory
    topo = "Nobel_EU"  # Replace with your topological name
    G, _, _ = read_graph(directory, topo)
    
    # Probability map (example, replace with actual data)
    proba_map = {node: 0.9 for node in G.nodes()}
    
    # Relabel the graph to a dictionary format
    G_relabel, proba_map_relabel, _ = relabel_graph_A_dict(G, proba_map)


    # Node pairs for evaluation
    node_pairs = list(combinations(G_relabel.nodes(), 2))

    mincutsets = read_mincutset(directory, topo)
    
    mincutsets = [
                [[i + 1 for i in num] for num in sublist] for sublist in mincutsets
            ]
    
    
    
    time_start = time.time()
    results = pyrbd_plusplus.mcs.eval_avail_topo(node_pairs, proba_map_relabel, mincutsets)
    print(f"Evaluation time: {time.time() - time_start:.2f} seconds") 

    # Evaluate minimal cuts
    time_start = time.time()
    results = pyrbd_plusplus.mcs.eval_avail_topo_parallel(node_pairs, proba_map_relabel, mincutsets)
    print(f"Evaluation time: {time.time() - time_start:.2f} seconds")