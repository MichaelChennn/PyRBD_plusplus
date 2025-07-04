import os
import pickle as pkl
import networkx as nx
import pandas as pd
import ast
from itertools import combinations
from tqdm import tqdm

import pyrbd_plusplus._core as pyrbd

from pyrbd_plusplus.algorithms import minimalcuts_optimized
from pyrbd_plusplus.algorithms import minimalpaths
from pyrbd_plusplus.utils import relabel_boolexpr_to_mathexpr

# Load data from a pickle file
def read_graph(directory, top):
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
        expression = build.rbd_bindings.MCSToProbaset(src, dst, mincutsets[i])

        # Append the boolean expression data
        boolean_expressions_data.append(
            {
                "source": src,
                "target": dst,
                "boolean_expression": relabel_boolexpr_to_mathexpr(expression),
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
        expression = build.rbd_bindings.pathSetToProbaset(src, dst, pathsets[i])

        # Append the boolean expression data
        boolean_expressions_data.append(
            {
                "source": src,
                "target": dst,
                "boolean_expression": relabel_boolexpr_to_mathexpr(expression),
                "length": len(expression),
            }
        )

    # Save boolean expressions to CSV
    boolean_expressions_df = pd.DataFrame(boolean_expressions_data)
    boolean_expressions_df.to_csv(
        os.path.join(directory, f"BoolExprPS_{top}.csv"), index=False
    )
    