import os
import sys
import time
import pandas as pd
from loguru import logger
from itertools import combinations
from tqdm import tqdm
from enum import Enum
import multiprocessing

from pyrbd_plusplus import *
from pyrbd_plusplus.algorithms.availability.pyrbd import (
    process_topology,
)
import pyrbd_plusplus._core.pyrbd_plusplus_cpp as cpp

# Define results directory
RESULTS_DIR = "benchmark/results"


# Define the Algorithm Enum for different algorithms
class Algorithm(Enum):
    PYRBD = "pyrbd"
    PYRBDPP_SDP = "pyrbdpp_sdp"
    PYRBDPP_PATHSET = "pyrbdpp_pathset"
    PYRBDPP_MCS = "pyrbdpp_mcs"
    PATHSET = "pathset"
    MCS = "mcs"
    OPTI_MCS = "opti_mcs"


def benchmark_mcs_pathset_per_flow(topology, skip_algorithms=[]):
    logger.info(
        f"===============Running Benchmark MCS and Pathset for {topology}==============="
    )

    # Load the topology
    logger.info(f"Loading topology: {topology}")
    G, _, _ = read_graph(f"topologies/{topology}", topology)

    # Get all node pairs
    node_pairs = list(combinations(list(G.nodes), 2))

    # Set the precision
    precision = 9
    logger.info(f"Setting the comparsion precision to {precision}")

    # Save the results each node pair
    result_per_flow = []

    # Save the results for each node pair
    mcs_time = 0
    mcs_opti_time = 0
    pathset_time = 0
    mcs_len = 0
    pathset_len = 0

    for src, dst in tqdm(
        node_pairs, desc=f"Running Benchmark MCS Pathset for {topology}"
    ):
        # Measure the MCS time taken
        if Algorithm.MCS not in skip_algorithms:
            mcs_time_start = time.time()
            minimalcuts(G, src, dst)
            mcs_time = time.time() - mcs_time_start

        # Measure the Optimized MCS time taken
        if Algorithm.OPTI_MCS not in skip_algorithms:
            mcs_opti_time_start = time.time()
            set_result = minimalcuts_optimized(G, src, dst)
            mcs_opti_time = time.time() - mcs_opti_time_start
            mcs_len = len(set_result)

        # Measure the Pathset time taken
        if Algorithm.PATHSET not in skip_algorithms:
            pathset_time_start = time.time()
            set_result = minimalpaths(G, src, dst)
            pathset_time = time.time() - pathset_time_start
            pathset_len = len(set_result)

        # Round the time taken
        mcs_time = round(mcs_time, precision)
        mcs_opti_time = round(mcs_opti_time, precision)
        pathset_time = round(pathset_time, precision)

        # Save the results for each node pair
        result_per_flow.append(
            {
                "Source": src,
                "Target": dst,
                "Pathset Time (Second)": pathset_time,
                "MCS Optimized Time (Second)": mcs_opti_time,
                "MCS Time (Second)": mcs_time,
                "Pathset Length": pathset_len,
                "MCS Length": mcs_len,
            }
        )
    logger.info(f"Benchmark MCS and Pathset for {topology} completed.")
    # Return the results as a DataFrame
    return pd.DataFrame(result_per_flow)


def benchmark_avail_per_flow(topology, skip_algorithms=[]):
    logger.info(
        f"===============Running Benchmark Availability Per Flow for {topology}==============="
    )

    # Warning for the user
    if topology in ["Nobel_EU"]:
        logger.warning(
            f"PyRBD is not scalable for topology: {topology}. The programm maybe be stucked."
        )
    elif topology in ["Germany_50"]:
        logger.warning(
            f"PyRBD, PyRBD++ MCS, PyRBD++ Pathset are not scalable for topology: {topology}. The programm maybe be stucked."
        )
        logger.warning(
            f"Please use the PyRBD++ SDP for this topology. The PyRBD++ SDP is scalable for this topology."
        )

    # Set the precision
    precision = 9
    logger.info(f"Setting the comparsion precision to {precision}")

    # Load the topology
    logger.info(f"Loading topology: {topology}")
    G, _, _ = read_graph(f"topologies/{topology}", topology)

    # Create the availability dictionaries
    A_dic = {i: 0.9 for i in range(G.number_of_nodes())}

    # Relabel the nodes of G and A_dic
    G_relabel, A_dic, relabel_mapping = relabel_graph_A_dict(G, A_dic)

    # Get all the node pairs
    node_pairs = list(combinations(list(G_relabel.nodes), 2))

    # Load the mincutsets
    try:
        logger.info("Loading mincutsets for " + topology)
        mincutsets, mincutsets_len = read_mincutset(f"topologies/{topology}", topology)
        mincutsets = [
            [[i + 1 for i in num] for num in sublist] for sublist in mincutsets
        ]
        mincutsets = [
            mincutset[2:] if len(mincutset) > 2 else mincutset
            for mincutset in mincutsets
        ]  # Filter out single node mincutsets
    except FileNotFoundError:
        logger.error(
            f"Mincutsets file for {topology} not found. Please generate the mincutsets first."
        )
        mincutsets = None

    # Load the pathsets
    try:
        logger.info("Loading pathsets for " + topology)
        pathsets, pathsets_len = read_pathset(f"topologies/{topology}", topology)
        pathsets = [[[i + 1 for i in num] for num in sublist] for sublist in pathsets]
    except FileNotFoundError:
        logger.error(
            f"Pathsets file for {topology} not found. Please generate the pathsets first."
        )
        pathsets = None

    # Save the results each node pair
    result_per_flow = []

    # Evaluate the availability with CPP
    for i, (src, dst) in tqdm(
        enumerate(node_pairs),
        total=len(node_pairs),
        desc=f"Running Benchmark Availability for {topology}",
    ):
        # Save the results per flow
        pyrbdpp_sdp_time = 0
        pyrbdpp_pathset_time = 0
        pyrbdpp_mcs_time = 0
        pyrbd_time = 0

        # Measure the PyRBD++ MCS and PyRBD time taken
        if mincutsets is not None:
            # Measure the PyRBD++ MCS time taken
            if Algorithm.PYRBDPP_MCS not in skip_algorithms:
                pyrbdpp_mcs_time_start = time.time()
                cpp.mcs.eval_avail(src, dst, A_dic, mincutsets[i])
                pyrbdpp_mcs_time = time.time() - pyrbdpp_mcs_time_start

            # Measure the PyRBD time taken
            if Algorithm.PYRBD not in skip_algorithms:
                time_start_pyrbd = time.time()
                process_topology(G_relabel, src, dst, A_dic, mincutsets=mincutsets[i])
                pyrbd_time = time.time() - time_start_pyrbd

        # Measure the PyRBD++ Pathset time and PyRBD++ SDP time taken
        if pathsets is not None:
            if Algorithm.PYRBDPP_PATHSET not in skip_algorithms:
                pyrbdpp_pathset_time_start = time.time()
                cpp.pathset.eval_avail(src, dst, A_dic, pathsets[i])
                pyrbdpp_pathset_time = time.time() - pyrbdpp_pathset_time_start

            # Measure the CPP SDP time taken
            if Algorithm.PYRBDPP_SDP not in skip_algorithms:
                pyrbdpp_sdp_time_start = time.time()
                cpp.sdp.eval_avail(src, dst, A_dic, pathsets[i])
                pyrbdpp_sdp_time = time.time() - pyrbdpp_sdp_time_start

        # Round the time taken to the specified precision
        pyrbdpp_sdp_time = round(pyrbdpp_sdp_time, precision)
        pyrbdpp_pathset_time = round(pyrbdpp_pathset_time, precision)
        pyrbdpp_mcs_time = round(pyrbdpp_mcs_time, precision)
        pyrbd_time = round(pyrbd_time, precision)

        # Append the results to the result_per_flow list
        result_per_flow.append(
            {
                "Source": src - 1,
                "Target": dst - 1,
                "PYRBD++ SDP Time (Second)": pyrbdpp_sdp_time,
                "PYRBD++ Pathset Time (Second)": pyrbdpp_pathset_time,
                "PYRBD++ MCS Time (Second)": pyrbdpp_mcs_time,
                "PYRBD Time (Second)": pyrbd_time,
                "Pathset Length": pathsets_len[i] if pathsets is not None else 0,
                "MCS Length": mincutsets_len[i] if mincutsets is not None else 0,
            }
        )

    logger.info(f"Benchmark Availability Per Flow for {topology} completed.")
    # Return the results as a DataFrame
    return pd.DataFrame(result_per_flow)


def benchmark_summary(topology, skip_algorithms=[]):
    logger.info(
        f"===============Running Benchmark Summary for {topology}==============="
    )
    # Load the topology
    logger.info(f"Loading topology: {topology}")
    G, _, _ = read_graph(f"topologies/{topology}", topology)

    # Read the avail per-flow results
    df_avail_per_flow = pd.read_csv(f"{RESULTS_DIR}/Avail_Per_Flow_{topology}.csv")
    
    # Read the mcs pathset per-flow results
    df_mcs_pathset_per_flow = pd.read_csv(
        f"{RESULTS_DIR}/MCS_Pathset_Per_Flow_{topology}.csv"
    )
    
    # Read three boolean expression files
    df_boolexpr_mcs = pd.read_csv(
        f"topologies/{topology}/BoolExprMCS_{topology}.csv"
    )
    df_boolexpr_pathset = pd.read_csv(
        f"topologies/{topology}/BoolExprPS_{topology}.csv"
    )
    df_boolexpr_sdp = pd.read_csv(
        f"topologies/{topology}/BoolExprSDP_{topology}.csv"
    )
    
    # Calculate the total time for PYRBD++ SDP, Pathset, MCS, and PYRBD
    pyrbdpp_sdp_time_total = df_avail_per_flow["PYRBD++ SDP Time (Second)"].sum()
    pyrbdpp_pathset_time_total = df_avail_per_flow[
        "PYRBD++ Pathset Time (Second)"
    ].sum()
    pyrbdpp_mcs_time_total = df_avail_per_flow["PYRBD++ MCS Time (Second)"].sum()
    pyrbd_time_total = df_avail_per_flow["PYRBD Time (Second)"].sum()

    # Calculate the total time for Pathset and MCS
    pathset_time_total = df_mcs_pathset_per_flow["Pathset Time (Second)"].sum()
    mcs_optimized_time_total = df_mcs_pathset_per_flow[
        "MCS Optimized Time (Second)"
    ].sum()
    mcs_time_total = df_mcs_pathset_per_flow["MCS Time (Second)"].sum()

    # Calculate the total time for some algorithms: PYRBD++ SDP + Pathset, PYRBD++ Pathset + Pathset, PYRBD++ MCS + Optimized MCS, PYRBD + MCS
    pyrbdpp_sdp_pathset_time_total = (
        0
        if pyrbdpp_sdp_time_total == 0.0 or pathset_time_total == 0.0
        else pyrbdpp_sdp_time_total + pathset_time_total
    )
    pyrbdpp_pathset_pathset_time_total = (
        0
        if pyrbdpp_pathset_time_total == 0.0 or pathset_time_total == 0.0
        else pyrbdpp_pathset_time_total + pathset_time_total
    )
    pyrbdpp_mcs_optimized_time_total = (
        0
        if pyrbdpp_mcs_time_total == 0.0 or mcs_optimized_time_total == 0.0
        else pyrbdpp_mcs_time_total + mcs_optimized_time_total
    )
    pyrbd_mcs_time_total = (
        0
        if pyrbd_time_total == 0.0 or mcs_time_total == 0.0
        else pyrbd_time_total + mcs_time_total
    )

    # Create the summary result dictionary
    summary_result = {
        "Topology": topology,
        "Nodes": len(G.nodes()),
        "Edges": len(G.edges()),
        "PYRBD++ SDP Time (Second)": pyrbdpp_sdp_time_total,
        "PYRBD++ Pathset Time (Second)": pyrbdpp_pathset_time_total,
        "PYRBD++ MCS Time (Second)": pyrbdpp_mcs_time_total,
        "PYRBD Time (Second)": pyrbd_time_total,
        "Pathset Time (Second)": pathset_time_total,
        "Optimized MCS Time (Second)": mcs_optimized_time_total,
        "Original MCS Time (Second)": mcs_time_total,
        "(Comb) PYRBD++ SDP Time (Second)": pyrbdpp_sdp_pathset_time_total,
        "(Comb) PYRBD++ Pathset Time (Second)": pyrbdpp_pathset_pathset_time_total,
        "(Comb) PYRBD++ MCS Time (Second)": pyrbdpp_mcs_optimized_time_total,
        "(Comb) PYRBD Time (Second)": pyrbd_mcs_time_total,
        "Pathset Length": df_avail_per_flow["Pathset Length"].sum(),
        "MCS Length": df_avail_per_flow["MCS Length"].sum(),
        "Bool Expr MCS Length": df_boolexpr_mcs["length"].sum(),
        "Bool Expr Pathset Length": df_boolexpr_pathset["length"].sum(),
        "Bool Expr SDP Length": df_boolexpr_sdp["length"].sum(),
    }

    logger.info(f"Benchmark Availability Summary for {topology} completed.")
    # Return the summary result as a DataFrame
    return pd.DataFrame([summary_result])


def benchmark_mcs_pathset_parallel(topology, skip_algorithms=[]):
    logger.info(
        f"===============Running Benchmark MCS and Pathset Parallel for {topology}==============="
    )

    # Load the topology
    logger.info(f"Loading topology: {topology}")
    G, _, _ = read_graph(f"topologies/{topology}", topology)

    # Get all node pairs
    node_pairs = list(combinations(list(G.nodes), 2))

    # Set the precision
    precision = 9
    logger.info(f"Setting the comparsion precision to {precision}")

    # Get the number of CPUs available
    num_cpus = multiprocessing.cpu_count()
    logger.info(f"Number of CPUs available: {num_cpus}")

    # Save the time taken for each algorithm
    mcs_time = 0
    mcs_opti_time = 0
    pathset_time = 0

    # Original MCS
    if Algorithm.MCS not in skip_algorithms:
        logger.info(f"Evaluating MCS Parallel")
        mcs_time_start = time.time()
        pool = multiprocessing.Pool(processes=num_cpus)
        pool.starmap(
            minimalcuts_optimized,
            [(G, pair[0], pair[1]) for pair in node_pairs],
        )
        pool.close()
        pool.join()
        mcs_time = time.time() - mcs_time_start

    # Optimized MCS
    if Algorithm.OPTI_MCS not in skip_algorithms:
        logger.info(f"Evaluating Optimized MCS Parallel")
        mcs_opti_time_start = time.time()
        pool = multiprocessing.Pool(processes=num_cpus)
        pool.starmap(
            minimalcuts_optimized,
            [(G, pair[0], pair[1]) for pair in node_pairs],
        )
        pool.close()
        pool.join()
        mcs_opti_time = time.time() - mcs_opti_time_start

    # Pathset
    if Algorithm.PATHSET not in skip_algorithms:
        logger.info(f"Evaluating Pathset Parallel")
        pathset_time_start = time.time()
        pool = multiprocessing.Pool(processes=num_cpus)
        pool.starmap(
            minimalpaths,
            [(G, pair[0], pair[1]) for pair in node_pairs],
        )
        pool.close()
        pool.join()
        pathset_time = time.time() - pathset_time_start

    # Round the time taken to the specified precision
    mcs_time = round(mcs_time, precision)
    mcs_opti_time = round(mcs_opti_time, precision)
    pathset_time = round(pathset_time, precision)

    # Save the results for each algorithm
    results = {
        "Topology": topology,
        "Pathset Parallel Time (Second)": pathset_time,
        "MCS Optimized Parallel Time (Second)": mcs_opti_time,
        "MCS Parallel Time (Second)": mcs_time,
    }
    logger.info(f"Benchmark MCS and Pathset Parallel for {topology} completed.")
    # Return the results as a DataFrame
    return pd.DataFrame([results])


def benchmark_avail_parallel(topology, skip_algorithms=[]):
    logger.info(
        f"===============Running Benchmark Availability Parallel for {topology}==============="
    )

    # Set the precision
    precision = 9
    logger.info(f"Setting the comparsion precision to {precision}")

    # Load the topology
    logger.info(f"Loading topology: {topology}")
    G, _, _ = read_graph(f"topologies/{topology}", topology)

    # Create the availability dictionaries
    A_dic = {i: 0.9 for i in range(G.number_of_nodes())}

    # Relabel the nodes of G and A_dic
    G_relabel, A_dic_relabel, relabel_mapping = relabel_graph_A_dict(G, A_dic)

    # Get all the node pairs
    node_pairs = list(combinations(list(G_relabel.nodes), 2))

    # Save the time taken for each algorithm
    pyrbdpp_sdp_time = 0
    pyrbdpp_pathset_time = 0
    pyrbdpp_mcs_time = 0
    pyrbd_time = 0

    # Save the lengths of the pathsets and mincutsets
    pathset_len = 0
    mcs_len = 0

    # Load the mincutsets
    try:
        logger.info(f"Loading mincutsets for {topology}")
        mincutsets, mincutsets_len = read_mincutset(f"topologies/{topology}", topology)
        mincutsets = [
            [[i + 1 for i in num] for num in sublist] for sublist in mincutsets
        ]
        mincutsets = [
            mincutset[2:] if len(mincutset) > 2 else mincutset
            for mincutset in mincutsets
        ]  # Filter out single node mincutsets
        mcs_len = sum(mincutsets_len)
    except FileNotFoundError:
        logger.error(
            f"Mincutsets file for {topology} not found. Please generate the mincutsets first. Check the pyrbd_plusplus.datasets module for preparing the mincutsets."
        )
        mincutsets = None

    # Load the pathsets
    try:
        logger.info(f"Loading pathsets for {topology}")
        pathsets, pathsets_len = read_pathset(f"topologies/{topology}", topology)
        pathsets = [[[i + 1 for i in num] for num in sublist] for sublist in pathsets]
        logger.debug(f"Loaded pathsets: {pathsets}")
        pathset_len = sum(pathsets_len)
    except FileNotFoundError:
        logger.error(
            f"Pathsets file for {topology} not found. Please generate the pathsets first. Check the pyrbd_plusplus.datasets module for preparing the pathsets."
        )
        pathsets = None

    # Get the number of CPUs available
    num_cpus = multiprocessing.cpu_count()
    logger.info(f"Number of CPUs available: {num_cpus}")

    # Evaluate the availability with PyRBD++ SDP
    if Algorithm.PYRBDPP_SDP not in skip_algorithms and pathsets is not None:
        logger.info(f"Evaluating availability PyRBD++ SDP Parallel")
        logger.debug(f"Availability dictionary: {A_dic_relabel}")
        pyrbdpp_sdp_time_start = time.time()
        cpp.sdp.eval_avail_topo_parallel(node_pairs, A_dic_relabel, pathsets)
        pyrbdpp_sdp_time = time.time() - pyrbdpp_sdp_time_start

    # Evaluate the availability with PyRBD++ Pathset
    if Algorithm.PYRBDPP_PATHSET not in skip_algorithms and pathsets is not None:
        logger.info(f"Evaluating availability PyRBD++ Pathset Parallel")
        pyrbdpp_pathset_time_start = time.time()
        cpp.pathset.eval_avail_topo_parallel(node_pairs, A_dic_relabel, pathsets)
        pyrbdpp_pathset_time = time.time() - pyrbdpp_pathset_time_start

    # Evaluate the availability with PyRBD++ MCS
    if Algorithm.PYRBDPP_MCS not in skip_algorithms and mincutsets is not None:
        logger.info(f"Evaluating availability PyRBD++ MCS Parallel")
        pyrbdpp_mcs_time_start = time.time()
        cpp.mcs.eval_avail_topo_parallel(node_pairs, A_dic_relabel, mincutsets)
        pyrbdpp_mcs_time = time.time() - pyrbdpp_mcs_time_start

    # Evaluate the availability with PyRBD
    if Algorithm.PYRBD not in skip_algorithms and mincutsets is not None:
        logger.info(f"Evaluating availability PyRBD Parallel")

        # Create a multiprocessing pool
        pyrbd_time_start = time.time()
        pool = multiprocessing.Pool(processes=num_cpus)
        pool.starmap(
            process_topology,
            [
                (G_relabel, pair[0], pair[1], A_dic_relabel, mincutset)
                for pair, mincutset in zip(node_pairs, mincutsets)
            ],
        )
        pool.close()
        pool.join()
        pyrbd_time = time.time() - pyrbd_time_start

    # Round the time taken to the specified precision
    pyrbdpp_sdp_time = round(pyrbdpp_sdp_time, precision)
    pyrbdpp_pathset_time = round(pyrbdpp_pathset_time, precision)
    pyrbdpp_mcs_time = round(pyrbdpp_mcs_time, precision)
    pyrbd_time = round(pyrbd_time, precision)

    # Create the summary result
    summary_result_parallel = {
        "Topology": topology,
        "PYRBD++ SDP Parallel Time (Second)": pyrbdpp_sdp_time,
        "PYRBD++ Pathset Parallel Time (Second)": pyrbdpp_pathset_time,
        "PYRBD++ MCS Parallel Time (Second)": pyrbdpp_mcs_time,
        "PYRBD Parallel Time (Second)": pyrbd_time,
        "Pathset Length": pathset_len,
        "MCS Length": mcs_len,
    }

    # Return the summary result as a DataFrame
    return pd.DataFrame([summary_result_parallel])


def benchmark_parallel_summary(topology, skip_algorithms=[]):
    logger.info(
        f"===============Running Benchmark Parallel Summary for {topology}==============="
    )

    # Load the topology
    logger.info(f"Loading topology: {topology}")
    G, _, _ = read_graph(f"topologies/{topology}", topology)

    # Get the number of CPUs available
    num_cpus = multiprocessing.cpu_count()
    logger.info(f"Number of CPUs available: {num_cpus}")
    
    # Read the avail parallel results
    df_avail_parallel = pd.read_csv(
        f"{RESULTS_DIR}/Avail_Parallel.csv"
    )
    
    # Read the mcs pathset parallel results
    df_mcs_pathset_parallel = pd.read_csv(
        f"{RESULTS_DIR}/MCS_Pathset_Parallel.csv"
    )
    
    # Handle the specific topology rows
    df_avail_parallel = df_avail_parallel[df_avail_parallel["Topology"] == topology]
    df_mcs_pathset_parallel = df_mcs_pathset_parallel[df_mcs_pathset_parallel["Topology"] == topology]
    
    # Extract times for PYRBD++ SDP, Pathset, MCS, and PYRBD
    pyrbdpp_sdp_parallel_time = df_avail_parallel["PYRBD++ SDP Parallel Time (Second)"].sum()
    pyrbdpp_pathset_parallel_time = df_avail_parallel["PYRBD++ Pathset Parallel Time (Second)"].sum()
    pyrbdpp_mcs_parallel_time = df_avail_parallel["PYRBD++ MCS Parallel Time (Second)"].sum()
    pyrbd_parallel_time = df_avail_parallel["PYRBD Parallel Time (Second)"].sum()

    # Extract times for Pathset, Optimized MCS and MCS
    pathset_parallel_time = df_mcs_pathset_parallel["Pathset Parallel Time (Second)"].sum()
    mcs_optimized_parallel_time = df_mcs_pathset_parallel["MCS Optimized Parallel Time (Second)"].sum()
    mcs_parallel_time = df_mcs_pathset_parallel["MCS Parallel Time (Second)"].sum()

    # Calculate combined times
    pyrbdpp_sdp_pathset_parallel_time = (
        0
        if pyrbdpp_sdp_parallel_time == 0.0 or pathset_parallel_time == 0.0
        else pyrbdpp_sdp_parallel_time + pathset_parallel_time
    )
    pyrbdpp_pathset_pathset_parallel_time = (
        0
        if pathset_parallel_time == 0.0 or pathset_parallel_time == 0.0
        else pathset_parallel_time + pathset_parallel_time
    )
    pyrbdpp_mcs_optimized_parallel_time = (
        0
        if mcs_parallel_time == 0.0 or mcs_optimized_parallel_time == 0.0
        else mcs_parallel_time + mcs_optimized_parallel_time
    )
    pyrbd_mcs_parallel_time = (
        0
        if pyrbd_parallel_time == 0.0 or mcs_parallel_time == 0.0
        else pyrbd_parallel_time + mcs_parallel_time
    )
    
    result = {
        "Topology": topology,
        "Nodes": len(G.nodes()),
        "Edges": len(G.edges()),
        "Number of CPUs": num_cpus,
        "PYRBD++ SDP Parallel Time (Second)": pyrbdpp_sdp_parallel_time,
        "PYRBD++ Pathset Parallel Time (Second)": pyrbdpp_pathset_parallel_time,
        "PYRBD++ MCS Optimized Parallel Time (Second)": pyrbdpp_mcs_optimized_parallel_time,
        "PYRBD++ MCS Parallel Time (Second)": pyrbdpp_mcs_parallel_time,
        "PYRBD Parallel Time (Second)": pyrbd_parallel_time,
        "Pathset Parallel Time (Second)": pathset_parallel_time,
        "Optimized MCS Parallel Time (Second)": mcs_optimized_parallel_time,
        "Original MCS Parallel Time (Second)": mcs_parallel_time,
        "(Comb) PYRBD++ SDP Parallel Time (Second)": pyrbdpp_sdp_pathset_parallel_time,
        "(Comb) PYRBD++ Pathset Parallel Time (Second)": pyrbdpp_pathset_pathset_parallel_time,
        "(Comb) PYRBD++ MCS Parallel Time (Second)": pyrbdpp_mcs_optimized_parallel_time,
        "(Comb) PYRBD Parallel Time (Second)": pyrbd_mcs_parallel_time,
    }

    logger.info(f"Benchmark Availability Parallel Summary for {topology} completed.")
    # Return the summary result as a DataFrame
    return pd.DataFrame([result])


# Define the benchmark type registry
type_REGISTRY = {
    "mcs_pathset_per_flow": {
        "columns": [
            "Source",
            "Target",
            "Pathset Time (Second)",
            "MCS Optimized Time (Second)",
            "MCS Time (Second)",
            "Pathset Length",
            "MCS Length",
        ],
        "run_func": benchmark_mcs_pathset_per_flow,
        "write_mode": "per_flow",  # Mode for writing results
    },
    "avail_per_flow": {
        "columns": [
            "Source",
            "Target",
            "PYRBD++ SDP Time (Second)",
            "PYRBD++ Pathset Time (Second)",
            "PYRBD++ MCS Time (Second)",
            "PYRBD Time (Second)",
            "Pathset Length",
            "MCS Length",
        ],
        "run_func": benchmark_avail_per_flow,  # Function to run the benchmark for this type
        "write_mode": "per_flow",  # Mode for writing results
    },
    "mcs_pathset_parallel": {
        "columns": [
            "Topology",
            "Pathset Parallel Time (Second)",
            "MCS Optimized Parallel Time (Second)",
            "MCS Parallel Time (Second)",
        ],
        "run_func": benchmark_mcs_pathset_parallel,
        "write_mode": "whole_topology",  # Mode for writing results
    },
    "avail_parallel": {
        "columns": [
            "Topology",
            "PYRBD++ SDP Parallel Time (Second)",
            "PYRBD++ Pathset Parallel Time (Second)",
            "PYRBD++ MCS Parallel Time (Second)",
            "PYRBD Parallel Time (Second)",
            "Pathset Length",
            "MCS Length",
        ],
        "run_func": benchmark_avail_parallel,
        "write_mode": "whole_topology",  # Mode for writing results
    },
    "summary": {
        "columns": [
            "Topology",
            "Nodes",
            "Edges",
            "PYRBD++ SDP Time (Second)",
            "PYRBD++ Pathset Time (Second)",
            "PYRBD++ MCS Time (Second)",
            "PYRBD Time (Second)",
            "Pathset Time (Second)",
            "Optimized MCS Time (Second)",
            "Original MCS Time (Second)",
            "(Comb) PYRBD++ SDP Time (Second)",
            "(Comb) PYRBD++ Pathset Time (Second)",
            "(Comb) PYRBD++ MCS Time (Second)",
            "(Comb) PYRBD Time (Second)",
            "Pathset Length",
            "MCS Length",
            "Bool Expr MCS Length",
            "Bool Expr Pathset Length",
            "Bool Expr SDP Length", 
        ],
        "run_func": benchmark_summary,
        "write_mode": "whole_topology",  # Mode for writing results
    },
    "summary_parallel": {
        "columns": [
            "Topology",
            "Nodes",
            "Edges",
            "Number of CPUs",
            "PYRBD++ SDP Parallel Time (Second)",
            "PYRBD++ Pathset Parallel Time (Second)",
            "PYRBD++ MCS Parallel Time (Second)",
            "PYRBD Parallel Time (Second)",
            "Pathset Parallel Time (Second)",
            "Optimized MCS Parallel Time (Second)",
            "Original MCS Parallel Time (Second)",
            "(Comb) PYRBD++ SDP Parallel Time (Second)",
            "(Comb) PYRBD++ Pathset Parallel Time (Second)",
            "(Comb) PYRBD++ MCS Parallel Time (Second)",
            "(Comb) PYRBD Parallel Time (Second)",
        ],
        "run_func": benchmark_parallel_summary,
        "write_mode": "whole_topology",  # Mode for writing results
    },
}

def benchmark_scalability():
    """ This benchmark is only for Germany
    """

class BenchmarkTask:
    def __init__(self, topology, type, result_file, skip_algorithms=[]):
        self.topology = topology
        self.type = type
        self.result_file = result_file
        self.columns = type_REGISTRY[type]["columns"]
        self.run_func = type_REGISTRY[type]["run_func"]
        self.write_mode = type_REGISTRY[type]["write_mode"]
        self.skip_algorithms = skip_algorithms
        dirname = os.path.dirname(result_file)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)

    def run_and_save(self, overwrite=False, *args, **kwargs):
        result_df = self.run_func(
            self.topology, skip_algorithms=self.skip_algorithms, *args, **kwargs
        )

        # --- per-flow mode ---
        if self.write_mode == "per_flow":
            # For per-flow mode, save the entire DataFrame
            result_df.to_csv(self.result_file, index=False, float_format="%.9f")
            logger.info(
                f"{self.topology} {self.type} results saved to {self.result_file}"
            )
            return

        # --- whole topology mode ---
        # Check if file exists
        if os.path.exists(self.result_file):
            existing_df = pd.read_csv(self.result_file)

            # Check if topology already exists
            if not overwrite and self.topology in existing_df["Topology"].values:
                logger.info(
                    f"No overwriting: {self.topology} {self.type} results already exist in {self.result_file}, skipping."
                )
                return

            # Remove existing entry for this topology
            existing_df = existing_df[existing_df["Topology"] != self.topology]
        else:
            # Create empty DataFrame with correct columns
            existing_df = pd.DataFrame(columns=self.columns)

        # Ensure the result DataFrame has the correct columns
        for col in self.columns:
            if col not in result_df.columns:
                result_df[col] = None

        # Reorder columns to match self.columns
        result_df = result_df[self.columns]

        # Concatenate existing data with new data
        final_df = pd.concat([existing_df, result_df], ignore_index=True)

        # Save to file
        final_df.to_csv(self.result_file, index=False, float_format="%.9f")
        logger.info(
            f"{self.topology} {self.type} results saved to {self.result_file}"
        )


def benchmark_factory(tasks, overwrite=False):
    for task in tasks:
        task.run_and_save(overwrite=overwrite)


def register_tasks(mapping, result_dir=RESULTS_DIR):
    """
    mapping: dict, e.g. {"Abilene": ["avail_per_flow", "summary"], ...}
    result_dir: str, the directory where results will be saved.
    """
    tasks = []
    for topology, conf in mapping.items():
        types = conf.get("types", [])
        skip_algorithms = conf.get("skip_algos", {})
        for t in types:
            if isinstance(skip_algorithms, dict):
                skip_algo = skip_algorithms.get(t, [])
            else:
                skip_algo = skip_algorithms
            if t == "avail_per_flow":
                result_file = os.path.join(result_dir, f"Avail_Per_Flow_{topology}.csv")
                
            elif t == "avail_parallel":
                result_file = os.path.join(
                    result_dir, f"Avail_Parallel.csv"
                )
            elif t == "mcs_pathset_per_flow":
                result_file = os.path.join(
                    result_dir, f"MCS_Pathset_Per_Flow_{topology}.csv"
                )
            elif t == "mcs_pathset_parallel":
                result_file = os.path.join(
                    result_dir, f"MCS_Pathset_Parallel.csv"
                )
            elif t == "summary":
                result_file = os.path.join(result_dir, f"Summary_Series.csv")
                
            elif t == "summary_parallel":
                result_file = os.path.join(result_dir, f"Summary_Parallel.csv")
                
            else:  # For other types, use the topology name directly
                logger.warning(
                    f"Unknown benchmark type: {t}. Using topology name for result file."
                )

            tasks.append(
                BenchmarkTask(
                    topology=topology,
                    type=t,
                    result_file=result_file,
                    skip_algorithms=skip_algo,  # List of algorithms to skip
                )
            )
    return tasks

if __name__ == "__main__":
    # Example usage
    logger.remove()
    logger.add(sys.stderr, level="INFO", format="{message}", colorize=True)

    # Create task mapping
    task_mapping = {
        "Abilene": {
            "types": [
                "mcs_pathset_per_flow",
                "avail_per_flow",
                "summary",
                "mcs_pathset_parallel",
                "avail_parallel",
                "summary_parallel",
            ],
            "skip_algos": {},
        },
        "dfn-bwin": {
            "types": [
                "mcs_pathset_per_flow",
                "avail_per_flow",
                "summary",
                "mcs_pathset_parallel",
                "avail_parallel",
                "summary_parallel",
            ],
            "skip_algos": {},
        },
        "Germany_17": {
            "types": [
                "mcs_pathset_per_flow",
                "avail_per_flow",
                "summary",
                "mcs_pathset_parallel",
                "avail_parallel",
                "summary_parallel",
            ],
            "skip_algos": {},
        },
        "HiberniaUk": {
            "types": [
                "mcs_pathset_per_flow",
                "avail_per_flow",
                "summary",
                "mcs_pathset_parallel",
                "avail_parallel",
                "summary_parallel",
            ],
            "skip_algos": {},
        },
        "polska": {
            "types": [
                "mcs_pathset_per_flow",
                "avail_per_flow",
                "summary",
                "mcs_pathset_parallel",
                "avail_parallel",
                "summary_parallel",
            ],
            "skip_algos": {},
        },
        "Nobel_EU": {
            "types": [
                "mcs_pathset_per_flow",
                "avail_per_flow",
                "summary",
                "mcs_pathset_parallel",
                "avail_parallel",
                "summary_parallel",
            ],
            "skip_algos": {
                "avail_per_flow": [Algorithm.PYRBD],
                "mcs_pathset_per_flow": [Algorithm.MCS, Algorithm.OPTI_MCS],
                "avail_parallel": [Algorithm.PYRBD],
                "mcs_pathset_parallel": [Algorithm.MCS, Algorithm.OPTI_MCS],
            },
        },
    }

    # Register tasks
    tasks = register_tasks(task_mapping, result_dir=RESULTS_DIR)

    # Run the benchmark
    benchmark_factory(tasks, overwrite=True)
    logger.info("Benchmark completed successfully.")
