"""
PyRBD++ Library Demo
===================

This demo showcases all the main features of the PyRBD++ library for network reliability analysis.
The demo covers:
1. Basic usage with different algorithms (MCS, PathSet, SDP, PyRBD)
2. Single node-pair vs full topology evaluation
3. Sequential vs parallel processing
4. Small, medium, and large network examples
5. Performance comparisons and best practices

Author: 
License: MIT
"""
import time
from pyrbd_plusplus import read_graph, evaluate_availability

def print_section(title):
    print("\n" + "="*80)
    print(title)
    print("="*80 + "\n")

def print_subsection(title):
    print("\n" + "-"*80)
    print(title)
    print("-"*80)

def print_avail_result(result, title="Result", max_show=5):
    if isinstance(result, list):
        print(f"{title}: list of {len(result)} (src, dst, availability) tuples.")
        print("Sample results:")
        for item in result[:max_show]:
            print(f"  src={item[0]}, dst={item[1]}, availability={item[2]:.8f}")
        if len(result) > max_show:
            print("  ...")
    elif isinstance(result, tuple) and len(result) == 3:
        print(f"{title}: src={result[0]}, dst={result[1]}, availability={result[2]:.8f}")
    else:
        print(f"{title}: (unknown format) {result}")

print_section("Welcome to PyRBD++ - Fast Network Reliability Analysis Library!")

# Load a sample graph
print_subsection("Loading sample graph: Germany_17")
topo = "Germany_17"
G, _, _ = read_graph(f"topologies/{topo}", topo)
print(f"Graph {topo} loaded with {len(G.nodes())} nodes and {len(G.edges())} edges.")

# Load a sample graph another way
print_subsection("Alternative: Load graph by file path")
file_path = "topologies/Germany_17/Pickle_Germany_17.pickle"
G, _, _ = read_graph("", topo, file_path)

# Create a dictionary of node probabilities
print_subsection("Creating node probability dictionary (p=0.9)")
nodes_probabilities = {node: 0.9 for node in G.nodes()}
sample_prob = dict(list(nodes_probabilities.items())[:5])
print("Sample probabilities (first 5 nodes):", sample_prob)

# Evaluate availability for a single node pair
src_example = 0
dst_example = 1
print_section(f"Evaluate availability from node {src_example} to {dst_example} (Germany_17)")

for algo in ["mcs", "pathset", "sdp", "pyrbd"]:
    print_subsection(f"Algorithm: {algo.upper()}")
    start = time.time()
    avail = evaluate_availability(G, nodes_probabilities, src=src_example, dst=dst_example, algorithm=algo)
    t = time.time() - start
    print_avail_result(avail, f"Availability ({algo.upper()}, single pair)")
    print(f"Time: {t:.4f} s")

print_section("Evaluate full topology (all pairs) using SDP algorithm")
print("Note: Only show first 5 results, full result is a list of all node pairs.")
start = time.time()
avail_sdp_full = evaluate_availability(G, nodes_probabilities, algorithm="sdp")
t = time.time() - start
print_avail_result(avail_sdp_full, "Availability (SDP, all pairs)")
print(f"Time: {t:.4f} s")

print_subsection("SDP full topology, parallel processing")
start = time.time()
avail_sdp_full_parallel = evaluate_availability(G, nodes_probabilities, algorithm="sdp", parallel=True)
t = time.time() - start
print_avail_result(avail_sdp_full_parallel, "Availability (SDP, all pairs, parallel)")
print(f"Time: {t:.4f} s")

# Large graph demo
print_section("Large Topology Example: Germany_50")
large_topo = "Germany_50"
G_large, _, _ = read_graph(f"topologies/{large_topo}", large_topo)
nodes_probabilities_large = {node: 0.9 for node in G_large.nodes()}
print(f"Large graph loaded: {large_topo} ({len(G_large.nodes())} nodes, {len(G_large.edges())} edges)")

src_large = 0
dst_large = 1

print_subsection("SDP (single pair, large graph)")
start = time.time()
avail_large_sdp = evaluate_availability(G_large, nodes_probabilities_large, src=src_large, dst=dst_large, algorithm="sdp")
t = time.time() - start
print_avail_result(avail_large_sdp, "Availability (SDP, large graph, single pair)")
print(f"Time: {t:.4f} s")

print_subsection("SDP (single pair, large graph, parallel)")
start = time.time()
avail_large_sdp_parallel = evaluate_availability(G_large, nodes_probabilities_large, src=src_large, dst=dst_large, algorithm="sdp", parallel=True)
t = time.time() - start
print_avail_result(avail_large_sdp_parallel, "Availability (SDP, large graph, single pair, parallel)")
print(f"Time: {t:.4f} s")

print_section("Demo Completed!")
print("\nYou can customize algorithms, source/destination nodes, probability values, and graphs as needed.\n")

