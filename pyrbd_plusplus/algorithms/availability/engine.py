from itertools import combinations
from loguru import logger

from ...algorithms.sets import minimalcuts_optimized, minimalpaths
from ...utils import relabel_graph_A_dict, relabel_boolexpr_to_str, sdp_boolexpr_to_str
from .pyrbd import eval_avail_pyrbd, eval_avail_topo_pyrbd, eval_avail_topo_pyrbd_parallel
import pyrbd_plusplus._core.pyrbd_plusplus_cpp as cpp

# Algorithm Configuration
ALGORITHM_CONFIG = {
    'mcs': {
        'cpp_module': 'mcs',
        'problem_set_func': minimalcuts_optimized,
        'bool_expr_func': relabel_boolexpr_to_str,
        'to_set_func': 'to_probaset',
    },
    'pathset': {
        'cpp_module': 'pathset', 
        'problem_set_func': minimalpaths,
        'bool_expr_func': relabel_boolexpr_to_str,
        'to_set_func': 'to_probaset',
    },
    'sdp': {
        'cpp_module': 'sdp',
        'problem_set_func': minimalpaths,  # SDP is based on pathsets
        'bool_expr_func': sdp_boolexpr_to_str,
        'to_set_func': 'to_sdp_set',
    },
    'pyrbd': {
        'cpp_module': None,
        'problem_set_func': None,
        'bool_expr_func': None,
        'to_set_func': None,
    }
}

def to_boolean_expression(G, src, dst, algorithm='mcs'):
    """Convert a SUM of Disjoint Products (SDP) to a Boolean expression.
    
    Args:
        G (networkx.Graph): The graph representation of the topology.
        src (int): The source node index.
        dst (int): The destination node index.
        algorithm (str): The algorithm to use for generating the Boolean expression.
    
    Returns:
        bool_expr (str): A string representation of the Boolean expression.
    """
    if algorithm not in ALGORITHM_CONFIG:
        raise ValueError(f"Unsupported algorithm: {algorithm}. Choose from {list(ALGORITHM_CONFIG.keys())}.")
    
    if algorithm == 'pyrbd':
        raise NotImplementedError("Boolean expression generation for PyrBD is not implemented.")
    
    # Get algorithm configuration
    config = ALGORITHM_CONFIG[algorithm]
    cpp_module = getattr(cpp, config['cpp_module'])
    
    # Relabel
    G_relabel, _, relabel_mapping = relabel_graph_A_dict(G, {})
    src_relabel = relabel_mapping[src]
    dst_relabel = relabel_mapping[dst]
    
    # Get problem sets
    problem_sets = config['problem_set_func'](G_relabel, src_relabel, dst_relabel)
    
    # Convert to probability set or SDP set
    to_set_func = getattr(cpp_module, config['to_set_func'])
    result_set = to_set_func(src_relabel, dst_relabel, problem_sets)
    
    # Convert to boolean expression
    return config['bool_expr_func'](result_set)

def eval_single_pair(G, A_dict, src, dst, algorithm='mcs', parallel=False):
    """Evaluate the availability of a source and destination pair in the topology using MCS.
    
    Args:
        G (networkx.Graph): The graph representation.
        A_dict (dict): A dictionary mapping nodes to their availability.
        src (int): The source node index.
        dst (int): The destination node index.
        algorithm (str): The algorithm to use for evaluation ('mcs', 'pathset', or 'sdp').
        parallel (bool): Whether to use parallel evaluation if available.
    
    Raises:
        ValueError: If the specified algorithm does not support parallel evaluation.
    
    Returns:
        tuple: (src, dst, availability)
    """
    # Validate algorithm
    if algorithm not in ALGORITHM_CONFIG:
        raise ValueError(f"Unsupported algorithm: {algorithm}. Choose from {list(ALGORITHM_CONFIG.keys())}.")
    
    if algorithm == 'pyrbd':
        return eval_avail_pyrbd(G, src, dst, A_dict)
    
    # Get algorithm configuration
    config = ALGORITHM_CONFIG[algorithm]
    cpp_module = getattr(cpp, config['cpp_module'])
    
    # Relabel
    G_relabel, A_dict_relabeled, relabel_mapping = relabel_graph_A_dict(G, A_dict)
    src_relabel = relabel_mapping[src]
    dst_relabel = relabel_mapping[dst]
    
    # Get problem sets
    problem_sets = config['problem_set_func'](G_relabel, src_relabel, dst_relabel)
    
    # Choose evaluation method
    if parallel:
        if hasattr(cpp_module, 'eval_avail_parallel'):
            availability = cpp_module.eval_avail_parallel(
                src_relabel, dst_relabel, A_dict_relabeled, problem_sets
            )
        else:
            raise ValueError(f"Parallel evaluation not available for {algorithm} algorithm.")
    else:
        availability = cpp_module.eval_avail(
            src_relabel, dst_relabel, A_dict_relabeled, problem_sets
        )
    
    return (src, dst, availability)

def eval_topology(G, A_dict, algorithm='mcs', parallel=False):
    """Evaluate the availability for all pairs of nodes in the topology using MCS.
    
    Args:
        G (networkx.Graph): The graph representation.
        A_dict (dict): A dictionary mapping nodes to their availability.
        algorithm (str): The algorithm to use for evaluation ('mcs', 'pathset', or 'sdp').
        parallel (bool): Whether to use parallel evaluation if available.
        
    Raises:
        ValueError: If the specified algorithm does not support parallel evaluation.
    
    Returns:
        List[tuple]: A list of tuples, each containing (src, dst, availability).
    """
    # Validate algorithm
    if algorithm not in ALGORITHM_CONFIG:
        raise ValueError(f"Unsupported algorithm: {algorithm}. Choose from {list(ALGORITHM_CONFIG.keys())}.")
    
    if algorithm == 'pyrbd':
        return eval_avail_topo_pyrbd(G, A_dict) if not parallel else eval_avail_topo_pyrbd_parallel(G, A_dict)
    
    # Get algorithm configuration
    config = ALGORITHM_CONFIG[algorithm]
    cpp_module = getattr(cpp, config['cpp_module'])
    
    # Relabel
    G_relabel, A_dict_relabeled, relabel_mapping = relabel_graph_A_dict(G, A_dict)
    reverse_mapping = {v: k for k, v in relabel_mapping.items()}
    
    # Get all pairs
    node_pairs = list(combinations(G_relabel.nodes(), 2))
    
    # Get all problem sets
    problem_sets_list = [
        config['problem_set_func'](G_relabel, src, dst)
        for src, dst in node_pairs
    ]
    
    # Choose evaluation method
    if parallel:
        if hasattr(cpp_module, 'eval_avail_topo_parallel'):
            availability_lst = cpp_module.eval_avail_topo_parallel(
                node_pairs, A_dict_relabeled, problem_sets_list
            )
        else:
            raise ValueError(f"Parallel evaluation not available for {algorithm} algorithm.")
    else:
        availability_lst = cpp_module.eval_avail_topo(
            node_pairs, A_dict_relabeled, problem_sets_list
        )
    
    # Relabel results
    return [
        (reverse_mapping[src], reverse_mapping[dst], availability)
        for src, dst, availability in availability_lst
    ]
