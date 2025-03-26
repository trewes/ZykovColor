import os
import random
import math
import argparse
from pathlib import Path
from itertools import product
import networkx as nx


def generate_ER_graph(n, p):
    """Generate an Erdős–Rényi graph with n vertices and edge probability p."""
    assert (n > 0) and (0 < p < 1)
    edges = []
    for u in range(n):
        for v in range(u + 1, n):
            if random.random() <= p:
                edges.append((u, v))
    return edges


def generate_RGG(n, d):
    """Generate a random geometric graph U_{n, d}
        with n vertices randomly sampled from the unit disk and
        edges if their distance is less or equal to the threshold d.
        if d is of the form 1.x take the complement graph for d = 0.x"""
    assert (n > 0) and (0 < d < 2)
    complement = False
    if d > 1:
        d = round(d - 1, 15)
        complement = True
    points = [(random.uniform(0, 1), random.uniform(0, 1)) for _ in range(n)]
    edges = []
    for u in range(n):
        for v in range(u + 1, n):
            x1, y1 = points[u]
            x2, y2 = points[v]
            euclidean_distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
            if (not complement) and euclidean_distance <= d:
                edges.append((u, v))
            elif complement and euclidean_distance > d:
                edges.append((u, v))
    return edges


def generate_HGG(n, R=8, alpha=1.0, curvature=1.0):
    """ Provided by Jacob Miller
    Generates a hyperbolic geometric graph on n nodes (not guaranteed connected)
    n: Number of nodes in the graph
    R: the radius of the hyperbolic disk to sample from. R >> 1 to be distinct from Euclidean graph
    alpha: controls the degree distribution (large alpha -> longer tailed distribution)
    curvature: The curvature of the sampled space. Larger curvature -> fewer connections

    returns: list of edge tuples representing a graph

    Implemented according to https://arxiv.org/abs/1006.5169 and https://arxiv.org/abs/1503.05180.
    """
    from math import sinh, cosh, cos, acosh, pi
    import random
    assert alpha > 0
    assert curvature > 0

    quasi_uniform_factor = lambda u: acosh(1 + (cosh(alpha * R) - 1) * u) / alpha
    radii  = [quasi_uniform_factor(random.random()) for _ in range(n)]
    angles = [random.random() * (2*pi) for _ in range(n)]

    edges = []
    for u in range(n):
        for v in range(u + 1, n):
            r1,r2 = radii[u], radii[v]
            t1,t2 = angles[u], angles[v]
            r1 *= curvature
            r2 *= curvature

            hdist = cosh(r1) * cosh(r2) - sinh(r1) * sinh(r2) * cos(pi - abs(pi - abs(t1 - t2)))
            dist  = acosh(hdist) / curvature

            if dist <= R:
                edges.append((u,v))
    return edges


def generate_BA_graph(n, m):
    """Generate a Barabási–Albert graph with n vertices and m edges to attach from a new node to existing nodes."""
    assert (n > 0) and (m  < n)
    G = nx.barabasi_albert_graph(n, m)
    edges = list(G.edges())
    return edges



def write_graph_to_dimacs(filename, n, edges):
    """Write a graph in DIMACS format to a file."""
    with open(filename, "w") as file:
        file.write(f"c Random graph with {n} vertices\n")
        file.write(f"p edge {n} {len(edges)}\n")
        for u, v in edges:
            file.write(f"e {u + 1} {v + 1}\n")


def generate_graphs(output_dir, k, n, param, model):
    """Generate k graphs, write them in DIMACS format, and store them in a directory."""
    pname = str(param)
    if model == "ER": #ER graphs have the name ERn.p where 0.p is the probability parameter
        assert 0 < param < 1
        pname = str(param).split('.')[1] #get fractional part
    elif model == "RGG": #RGG graphs have the name RGGn.d where 0.d is the threshold parameter (*c for the complement)
        assert 0 < param < 2 and param != 1
        if param < 1:
            pname = str(param).split('.')[1]  # get fractional part
        elif param > 1:
            pname = str(round(param - 1, 15)).split('.')[1] + "c"

    dir_name = Path(output_dir) / f"{model}/{model}{n}.{pname}"
    if os.path.exists(dir_name):
        print(f"{dir_name} already exists. Skipping generation.")
        return
    Path(dir_name).mkdir(parents=True, exist_ok=True)
    for i in range(1, k + 1):
        if model == "ER":
            edges = generate_ER_graph(n, param)
        elif model == "RGG":
            edges = generate_RGG(n, param)
        elif model == "HGG":
            edges = generate_HGG(n, param)
        elif model == "BA":
            edges = generate_BA_graph(n, int(param))
        else:
            raise ValueError("Unsupported graph model: Use 'ER', 'RGG', 'HGG', or 'BA'")

        file_name = os.path.join(dir_name, f"{model}{n}.{pname}_inst{i}.col")
        write_graph_to_dimacs(file_name, n, edges)
    print(f"{k} Graphs with {n} vertices and parameter {pname} written to {dir_name}")


def generate_for_combinations(output_dir, k, n_set, param_set, model):
    """Meta-function to generate graphs for combinations of n and parameter."""
    for n, param in product(n_set, param_set):
        generate_graphs(output_dir, k, n, param, model)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate k random graphs for combinations of n and parameters, using specified models.",
        epilog="Example usage: python3 graph_gen.py -o ../Graphs/ --model ER -k 10 --n_set 100 125 150 --p_set 0.1 0.5 0.9"
    )
    parser.add_argument("-o", type=str, required=True, help="Parent directory where the graphs will be stored")
    parser.add_argument("-k", type=int, required=True, help="Number of graphs to generate for each combination")
    parser.add_argument("--n_set", type=int, nargs="+", required=True, help="List of vertex counts (n values)")
    parser.add_argument("--p_set", type=float, nargs="+", required=True,
                        help="List of parameters (edge probabilities for ER, threshold for RGG, radius for HGG, or preferential attachment for BA)")
    parser.add_argument("--model", type=str, required=True, choices=["ER", "RGG", "HGG", "BA"],
                        help="Graph model to use: 'ER' for Erdős–Rényi, 'RGG' for Random Geometric Graph, 'HGG' for Hyperbolic Geometric Graph, 'BA' for Barabási–Albert")

    args = parser.parse_args()
    print(args)
    generate_for_combinations(args.o, args.k, args.n_set, args.p_set, args.model)


# Example usage
# python graph_gen.py -o ../Graphs/random/ --model ER -k 100 --n_set 50 60 70 80 90 100 125 150 --p_set 0.1 0.3 0.5 0.7 0.9

