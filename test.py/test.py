from pyrbd_plusplus import read_graph, evaluate_availability
from itertools import combinations
import pandas as pd

def main():
    G, _, _ = read_graph("topologies/Germany_17", "Germany_17")
    A_dict = {node: 0.99 for node in G.nodes()}

    for src, dst in combinations(G.nodes(), 2):
        avail = evaluate_availability(G, A_dict, src=src, dst=dst, algorithm="sdp")
        print(f"Availability from {src} to {dst}: {avail}")

if __name__ == "__main__":
    main()
    topologies = ["Abilene", "dfn-bwin", "Germany_17", "HiberniaUk", "Nobel_EU", "polska"]
    for topology in topologies:
        avail_csv = f"topologies/{topology}/Availability_{topology}.csv"
        df = pd.read_csv(avail_csv)
        df["source"] = df["source"] - 1
        df["target"] = df["target"] - 1
        
        df.to_csv(avail_csv, index=False)
        print(f"Processed {avail_csv} successfully.")