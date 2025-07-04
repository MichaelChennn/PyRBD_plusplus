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
def relabel_boolexpr_to_mathexpr(bool_expr):
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