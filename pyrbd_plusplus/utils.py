import networkx as nx

# Relabel the nodes of G and A_dic
def relabel_graph_A_dict(G, A_dic):
    # Get the nodes of G
    nodes = list(G.nodes())
    # Create a mapping of the nodes to new labels
    relabel_mapping = {nodes[i]: i + 1 for i in range(len(nodes))}
    # Relabel the nodes of G
    G_relabel = nx.relabel_nodes(G, relabel_mapping)
    # Relabel the nodes in A_dict
    A_dic = {relabel_mapping[node]: value for node, value in A_dic.items()}
    return G_relabel, A_dic, relabel_mapping

# relabel boolean expression back and convert the boolean expression to a mathematical expression
def relabel_boolexpr_to_str(bool_expr):
    expr = "["
    bool_expr_len = len(bool_expr)
    for num_list in bool_expr:
        expr += "["
        list_len = len(num_list)
        for num in num_list:
            if num == -1:
                expr += "0"
            elif num > 0:
                expr += str(num - 1)
            elif num < 0:
                expr += str(num + 1)
            list_len -= 1
            if list_len > 0:
                expr += " * "
        expr += "]"
        bool_expr_len -= 1
        if bool_expr_len > 0:
            expr += " + "
    expr += "]"
    return expr

def sdp_boolexpr_to_str(sdp_lst_lst):
    sdp_str = ""
    sdp_lst_lst_len = len(sdp_lst_lst)
    for sdp_lst in sdp_lst_lst:
        sdp_lst_len = len(sdp_lst)
        sdp_str += "["
        for sdp in sdp_lst:

            if sdp.isComplementary():
                sdp_str += "-["
            else:
                sdp_str += "["

            sdp_elem_len = len(sdp.getSet())
            for sdp_elem in sdp.getSet():
                if sdp_elem_len == 1:
                    sdp_str += f"{sdp_elem}"
                else:
                    sdp_str += f"{sdp_elem} * "
                sdp_elem_len -= 1
            sdp_str += "]"

            sdp_lst_len -= 1

            if sdp_lst_len > 0:
                sdp_str += " * "

        sdp_str += "]"
        sdp_lst_lst_len -= 1
        if sdp_lst_lst_len > 0:
            sdp_str += " + "

    return sdp_str.strip()


def sdp_boolexpr_length(sdp_lst_lst):
    expr_len = 0
  
    for sdp_lst in sdp_lst_lst:
        
        sdp_lst_len = len(sdp_lst)
        
        expr_len += sdp_lst_len

    return expr_len
