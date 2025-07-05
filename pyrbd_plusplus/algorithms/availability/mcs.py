from .engine import to_boolean_expression, eval_single_pair, eval_topology

def mcs_to_boolexpr(G, src, dst):
    """Convert a Minimal Cut Set (MCS) to a Boolean expression."""
    return to_boolean_expression(G, src, dst, 'mcs')

def eval_avail_mcs(G, A_dict, src, dst):
    """Evaluate the availability of a Minimal Cut Set (MCS) for a single source-destination pair."""
    return eval_single_pair(G, A_dict, src, dst, 'mcs')

def eval_avail_topo_mcs(G, A_dict):
    """Evaluate the availability of a Minimal Cut Set (MCS) for the entire topology."""
    return eval_topology(G, A_dict, 'mcs')

def eval_avail_topo_mcs_parallel(G, A_dict):
    """Evaluate the availability of a Minimal Cut Set (MCS) for the entire topology in parallel."""
    return eval_topology(G, A_dict, 'mcs', parallel=True)