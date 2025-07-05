import os
import pickle as pkl
import pandas as pd
import ast
from itertools import combinations
from tqdm import tqdm

from .algorithms.sets import minimalcuts_optimized, minimalpaths
from .utils import relabel_boolexpr_to_str, sdp_boolexpr_to_str, sdp_boolexpr_length

import pyrbd_plusplus._core.pyrbd_plusplus_cpp as cpp

# Load data from a pickle file
def read_graph(directory, top, file_path=None):
    if file_path:
        with open(file_path, "rb") as handle:
            f = pkl.load(handle)
    else:
        with open(os.path.join(directory, "Pickle_" + top + ".pickle"), "rb") as handle:
            f = pkl.load(handle)
    G = f[0]
    pos = f[1]
    lable = f[2]

    return G, pos, lable

# Read mincutset from a csv file
def read_mincutset(directory, top):
    # Read the mincutset from a csv file
    df = pd.read_csv(os.path.join(directory, "Mincutset_" + top + ".csv"))

    # Convert the 'min-cutsets' column from string representation to actual list
    df["mincutsets"] = df["mincutsets"].apply(ast.literal_eval)

    # Extract the 'min-cutsets' column as a list
    mincutsets = df["mincutsets"].values.tolist()

    return mincutsets


# Read pathset from a csv file
def read_pathset(directory, top):
    # Read the pathset from a csv file
    df = pd.read_csv(os.path.join(directory, "Pathset_" + top + ".csv"))

    # Convert the 'pathsets' column from string representation to actual list
    df["pathsets"] = df["pathsets"].apply(ast.literal_eval)

    # Extract the 'pathsets' column as a list
    pathsets = df["pathsets"].values.tolist()

    return pathsets

# Evaluate the mincutset and save it to a csv file
def save_mincutset(directory, top):
    # Read the graph from the pickle file
    G, _, _ = read_graph(directory, top)

    # Get all node pairs
    node_pairs = list(combinations(G.nodes(), 2))

    # Initialize a list to store mincutsets data
    mincutsets_data = []

    # Iterate through each pair of nodes and find the minimal cut sets
    for src, dst in tqdm(node_pairs, desc=f"Saving Mincutsets for {top}", leave=False):
        mincutset = minimalcuts_optimized(G, src, dst)
        if mincutset:
            mincutsets_data.append(
                {
                    "source": src,
                    "target": dst,
                    "mincutsets": repr(mincutset),
                    "length": len(mincutset),
                }
            )

    # Save mincutsets to CSV
    mincutsets_df = pd.DataFrame(mincutsets_data)
    mincutsets_df.to_csv(os.path.join(directory, f"Mincutset_{top}.csv"), index=False)


# Evaluate the pathset and save it to a csv file
def save_pathset(directory, top):
    # Read the graph from the pickle file
    G, _, _ = read_graph(directory, top)

    # Get all node pairs
    node_pairs = list(combinations(G.nodes(), 2))

    # Initialize a list to store pathsets data
    pathsets_data = []

    # Iterate through each pair of nodes and find the simple paths
    for src, dst in tqdm(node_pairs, desc=f"Saving Pathsets for {top}", leave=False):
        pathset = minimalpaths(G, src, dst)
        if pathset:
            # Sort paths by length and for the same length, sort by lexicographical order
            pathset = sorted(pathset, key=lambda x: (len(x), x))
            pathsets_data.append(
                {
                    "source": src,
                    "target": dst,
                    "pathsets": repr(pathset),
                    "length": len(pathset),
                }
            )

    # Save pathsets to CSV
    pathsets_df = pd.DataFrame(pathsets_data)
    pathsets_df.to_csv(os.path.join(directory, f"Pathset_{top}.csv"), index=False)
    
    # Evaluate the boolean expression from the mincutsets and save it to a csv file
def save_boolean_expression_from_mincutset(directory, top):
    # Read the mincutset from the csv file
    mincutsets = read_mincutset(directory, top)

    # Read the graph from the pickle file
    G, _, _ = read_graph(directory, top)

    # Create node pairs
    node_pairs = list(combinations(G.nodes(), 2))

    # Initialize a list to store boolean expressions data
    boolean_expressions_data = []

    # Iterate through each mincutset and evaluate the boolean expression
    for i, (src, dst) in tqdm(
        enumerate(node_pairs), desc=f"Saving Boolean Expressions from MCS for {top}", 
        leave=False, total=len(node_pairs)
    ):
        # Relabel mincutset by adding 1 to ensure no 0 in the set
        mincutsets[i] = [[node + 1 for node in cutset] for cutset in mincutsets[i]]

        # Build the boolean expression from the mincutset
        expression = cpp.mcs.to_probaset(src, dst, mincutsets[i])

        # Append the boolean expression data
        boolean_expressions_data.append(
            {
                "source": src,
                "target": dst,
                "boolean_expression": relabel_boolexpr_to_str(expression),
                "length": len(expression),
            }
        )

    # Save boolean expressions to CSV
    boolean_expressions_df = pd.DataFrame(boolean_expressions_data)
    boolean_expressions_df.to_csv(
        os.path.join(directory, f"BoolExprMCS_{top}.csv"), index=False
    )


# Evaluate the boolean expression from the pathset and save it to a csv file
def save_boolean_expression_from_pathset(directory, top):
    # Read the pathset from the csv file
    pathsets = read_pathset(directory, top)

    # Read the graph from the pickle file
    G, _, _ = read_graph(directory, top)

    # Create node pairs
    node_pairs = list(combinations(G.nodes(), 2))

    # Initialize a list to store boolean expressions data
    boolean_expressions_data = []

    # Iterate through each pathset and evaluate the boolean expression
    for i, (src, dst) in tqdm(
        enumerate(node_pairs), desc=f"Saving Boolean Expressions from Pathset for {top}",
        leave=False, total=len(node_pairs)
    ):
        # Relabel pathset by adding 1 to ensure no 0 in the set
        pathsets[i] = [[node + 1 for node in path] for path in pathsets[i]]

        # Build the boolean expression from the pathset
        expression = cpp.pathset.to_probaset(src, dst, pathsets[i])
        
        # Append the boolean expression data
        boolean_expressions_data.append(
            {
                "source": src,
                "target": dst,
                "boolean_expression": relabel_boolexpr_to_str(expression),
                "length": len(expression),
            }
        )

    # Save boolean expressions to CSV
    boolean_expressions_df = pd.DataFrame(boolean_expressions_data)
    boolean_expressions_df.to_csv(
        os.path.join(directory, f"BoolExprPS_{top}.csv"), index=False
    )
    
def save_boolean_expression_from_sdp(directory, top):
    # Read the mincutset from the csv file
    pathsets = read_pathset(directory, top)

    # Read the graph from the pickle file
    G, _, _ = read_graph(directory, top)

    # Create node pairs
    node_pairs = list(combinations(G.nodes(), 2))

    # Initialize a list to store boolean expressions data
    boolean_expressions_data = []

    for i, (src, dst) in tqdm(
        enumerate(node_pairs), desc=f"Saving Boolean Expressions from SDP for {top}",
        leave=False, total=len(node_pairs)
    ):
        # Relabel pathset by adding 1 to ensure no 0 in the set
        pathsets[i] = [[node + 1 for node in path] for path in pathsets[i]]

        # Build the SDP set from the pathset
        expression = cpp.sdp.to_sdp_set(src, dst, pathsets[i])

        # Convert the SDP set to a string representation and calculate the expression length
        bool_expr_str = sdp_boolexpr_to_str(expression)

        # Append the result to the list
        boolean_expressions_data.append({
            "src": src,
            "dst": dst,
            "boolean_expression": bool_expr_str,
            "length": sdp_boolexpr_length(expression)
        })
    
    # Save the boolean expressions to a CSV file
    boolean_expressions_df = pd.DataFrame(boolean_expressions_data)
    boolean_expressions_df.to_csv(
        os.path.join(directory, f"BoolExprSDP_{top}.csv"), index=False
    )

def dataset_preparation(directory, topologies):
    """
    Prepare the dataset by saving pathsets, mincutsets, and boolean expressions.
    Args:
        directory (str): The directory where the datasets will be saved.
        topologies (list): A list of topology names to process.
    """
    for top in topologies:
        directory = f"topologies/{top}"
        save_pathset(directory, top)
        if top not in ["Germany_50", "Nobel_EU"]:
            save_mincutset(directory, top)
        if top not in ["Germany_50"]:
            save_boolean_expression_from_mincutset(directory, top)
            save_boolean_expression_from_pathset(directory, top)
            save_boolean_expression_from_sdp(directory, top)

if __name__ == "__main__":
    # Data set
    topologies = [
        "Abilene",
        "dfn-bwin",
        "Germany_17",
        "HiberniaUk",
        "polska",
        "Nobel_EU",
        "Germany_50",
    ]

    # Prepare the dataset
    dataset_preparation("topologies", topologies)
    