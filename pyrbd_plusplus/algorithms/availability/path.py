from .engine import to_boolean_expression, eval_single_pair, eval_topology

def pathset_to_boolexpr(G, src, dst):
    """Convert a pathset to a Boolean expression."""
    return to_boolean_expression(G, src, dst, 'pathset')

def eval_avail_pathset(G, A_dict, src, dst):
    """Evaluate the availability of a pathset for a single source-destination pair."""
    return eval_single_pair(G, A_dict, src, dst, 'pathset')

def eval_avail_topo_pathset(G, A_dict):
    """Evaluate the availability of a pathset for the entire topology."""
    return eval_topology(G, A_dict, 'pathset')

def eval_avail_topo_pathset_parallel(G, A_dict):
    """Evaluate the availability of a pathset for the entire topology in parallel."""
    return eval_topology(G, A_dict, 'pathset', parallel=True)
