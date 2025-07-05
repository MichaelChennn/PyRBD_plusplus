from .datasets import read_graph
from .algorithms.availability import eval_single_pair, eval_topology

def evaluate_availability(
    graph_or_filepath,
    nodes_probabilities,
    src=None,
    dst=None,
    algorithm="sdp",
    parallel=False,
):
    """Evaluate the availability of a network based on the given parameters.

    Args:
        graph_or_filepath (Union[str, nx.Graph]): Either a file path to the graph pickle data or a NetworkX graph object.
        nodes_probabilities (dict): A dictionary mapping node IDs to their availability probabilities.
        src (int, optional): The source node ID. Defaults to None.
        dst (int, optional): The destination node ID. Defaults to None.
        algorithm (str, optional): The algorithm to use for evaluation. Defaults to 'sdp'.
        parallel (bool, optional): Whether to evaluate in parallel. Defaults to False.
    """
    # Read the graph from the specified directory and topology
    if isinstance(graph_or_filepath, str):
        G, _, _ = read_graph("", "", graph_or_filepath)
    elif hasattr(graph_or_filepath, "nodes") and hasattr(graph_or_filepath, "edges"):
        # Assume graph_or_filepath is a networkx.Graph object
        G = graph_or_filepath
    else:
        raise TypeError(
            "Invalid graph input. Provide a file path or a NetworkX graph object."
        )
    # Validate the nodes_probabilities
    if not isinstance(nodes_probabilities, dict):
        raise TypeError(
            "nodes_probabilities must be a dictionary mapping node IDs to probabilities."
        )
    if len(nodes_probabilities) != len(G.nodes()):
        raise ValueError(
            "nodes_probabilities must contain probabilities for all nodes in the graph."
        )

    # Validate source and destination nodes
    if src is not None and dst is not None:
        print(
            f"Evaluating availability from node {src} to node {dst} using {algorithm} algorithm."
        )
        if src not in G.nodes() or dst not in G.nodes():
            raise ValueError(
                f"Source node {src} or destination node {dst} not found in the graph."
            )
        return eval_single_pair(G, nodes_probabilities, src, dst, algorithm, parallel)
    elif src is None and dst is None:
        print("Evaluating availability for all nodes")
        return eval_topology(G, nodes_probabilities, algorithm, parallel)
    else:
        raise ValueError(
            "Both source and destination nodes must be specified or neither."
        )
