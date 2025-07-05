from engine import to_boolean_expression, eval_single_pair, eval_topology

def sdp_to_boolexpr(G, src, dst):
    """Convert a Sum of Disjoint Products (SDP) to a Boolean expression."""
    return to_boolean_expression(G, src, dst, 'sdp')

def eval_avail_sdp(G, A_dict, src, dst):
    """Evaluate the availability of a Sum of Disjoint Products (SDP) for a single source-destination pair."""
    return eval_single_pair(G, A_dict, src, dst, 'sdp')

def eval_avail_topo_sdp(G, A_dict):
    """Evaluate the availability of a Sum of Disjoint Products (SDP) for the entire topology."""
    return eval_topology(G, A_dict, 'sdp')

def eval_avail_sdp_parallel(G, A_dict, src, dst):
    """Evaluate the availability of a Sum of Disjoint Products (SDP) for a single source-destination pair in parallel."""
    return eval_single_pair(G, A_dict, src, dst, 'sdp', parallel=True)

def eval_avail_topo_sdp_parallel(G, A_dict):
    """Evaluate the availability of a Sum of Disjoint Products (SDP) for the entire topology in parallel."""
    return eval_topology(G, A_dict, 'sdp', parallel=True)
