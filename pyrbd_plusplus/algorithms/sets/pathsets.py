import networkx as nx

# Function to get nodes that should not be visited again
def get_dontgo(G, visited):
    dontgo = []
    if len(visited) > 1:  # Dont care if only source in the visited list
        for index, node in enumerate(visited[:-1]):
            dontgo = dontgo + list(
                nx.neighbors(G, node)
            )  # Get all neighbors of previous node
    return dontgo


# USe this. removes unnecessary paths also
def all_simple_paths(G, source, target):
    visited = [source]
    stack = [iter(G[source])]
    while stack:
        children = stack[-1]
        child = next(children, None)
        dontgo = get_dontgo(G, visited)
        if child is None:
            stack.pop()
            visited.pop()
        elif child not in dontgo:
            if child == target:
                yield visited + [target]
            elif child not in visited:
                visited.append(child)
                stack.append(iter(G[child]))
                
# Function to get all minimal paths between source and destination
def minimalpaths(G, src, dst):
    return list(all_simple_paths(G, src, dst))