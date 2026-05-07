import os
import dendropy
import glob
import numpy as np

from snakemake_helpers import get_tool_params
from viz.utils import PROJECT_ROOT, load_snakemake_config_yaml, add_to_ordered_set, write_table

def tree_height(tree):
    """Sum of all branch lengths (tree height)."""
    s = 0.0
    root = tree.seed_node
    for node in tree.nodes():
        if node.edge and node.edge.length is not None and node is not root:
            s += node.edge.length
    return s

def num_leaves(tree: dendropy.Tree) -> int:
    """Number of taxa in the tree."""
    return len(tree.leaf_nodes())

def colless_index(tree: dendropy.Tree) -> float:
    """Colless index - normalized tree balance measure (0 = perfectly balanced, 1 = most unbalanced)."""
    n = len(tree.leaf_nodes())
    if n <= 2:
        return 0.0
    
    total = 0
    for node in tree.nodes():
        if not node.is_leaf() and len(node.child_nodes()) == 2:
            left_size = len(node.child_nodes()[0].leaf_nodes())
            right_size = len(node.child_nodes()[1].leaf_nodes())
            total += abs(left_size - right_size)
    
    max_colless = (n * (n - 1) * (n - 2)) / 4
    return total / max_colless if max_colless > 0 else 0.0

def leaf_distance_to_root(tree: dendropy.Tree):
    """Distance from all leaves to the root node (in branch length units)."""
    leaf_nodes = tree.leaf_nodes()
    if not leaf_nodes:
        return 0.0
    root = tree.seed_node
    total_distances = []
    for leaf in leaf_nodes:
        node = leaf
        while node is not root:
            total_distances.append(node.edge.length)
            node = node.parent_node
    d = {}
    d["avg_leaf_dist_to_root"] = sum(total_distances) / len(leaf_nodes)
    d["var_leaf_dist_to_root"] = np.var(total_distances)
    return d

def max_leaf_distance_to_root(tree: dendropy.Tree):
    """Maximum distance from any leaf to the root node (in branch length units)."""
    leaf_nodes = tree.leaf_nodes()
    if not leaf_nodes:
        return 0.0
    root = tree.seed_node
    max_distance = 0.0
    for leaf in leaf_nodes:
        node = leaf
        distance = 0.0
        while node is not root:
            distance += node.edge.length
            node = node.parent_node
        max_distance = max(max_distance, distance)
    
    return max_distance

def leaf_steps_to_root(tree: dendropy.Tree):
    """Average number of steps (edges) from all leaves to the root node."""
    root = tree.seed_node

    leaf_nodes = tree.leaf_nodes()
    if not leaf_nodes:
        return 0.0
    
    total_steps = []
    for leaf in leaf_nodes:
        steps = 0
        node = leaf
        while node is not root:
            steps += 1
            node = node.parent_node
        total_steps.append(steps)
    
    d = {}
    d["avg_leaf_steps_to_root"] = sum(total_steps) / len(leaf_nodes)
    d["var_leaf_steps_to_root"] = np.var(total_steps)
    return d

def max_leaf_steps_to_root(tree: dendropy.Tree) -> float:
    """Maximum number of steps (edges) from any leaf to the root node."""
    root = tree.seed_node

    leaf_nodes = tree.leaf_nodes()
    if not leaf_nodes:
        return 0.0
    
    max_steps = 0
    for leaf in leaf_nodes:
        steps = 0
        node = leaf
        while node is not root:
            steps += 1
            node = node.parent_node
        max_steps = max(max_steps, steps)
    
    return max_steps

def pairwise_leaf_step_distances(tree):
    """Compute pairwise distances (number of edges) between all pairs of leaves."""
    leaf_nodes = tree.leaf_nodes()
    if len(leaf_nodes) < 2:
        return None

    # Precompute the path from each leaf to the root
    leaf_to_root_paths = {}
    for leaf in leaf_nodes:
        path = []
        node = leaf
        while node is not None:
            path.append(node)
            node = node.parent_node
        leaf_to_root_paths[leaf] = path

    pairwise_distances = []
    for i in range(len(leaf_nodes)):
        for j in range(i + 1, len(leaf_nodes)):
            leaf1, leaf2 = leaf_nodes[i], leaf_nodes[j]
            path1, path2 = leaf_to_root_paths[leaf1], leaf_to_root_paths[leaf2]

            # Find the common ancestor
            # Start from the end of the paths (root) and move backwards until they diverge
            common_ancestor_index = 0
            while (common_ancestor_index < len(path1) and 
                   common_ancestor_index < len(path2) and 
                   path1[-(common_ancestor_index + 1)] == path2[-(common_ancestor_index + 1)]):
                common_ancestor_index += 1

            # Calculate distance as steps from each leaf to the common ancestor
            distance = (len(path1) - common_ancestor_index) + (len(path2) - common_ancestor_index)
            pairwise_distances.append(distance)
    d = {}
    d["avg_pairwise_leaf_step_distance"] = np.mean(pairwise_distances)
    d["var_pairwise_leaf_step_distance"] = np.var(pairwise_distances)
    return d

def main():
    config = load_snakemake_config_yaml()
    base_dir = os.path.join(PROJECT_ROOT, "results/sim/tree")
    os.makedirs(base_dir, exist_ok=True)
    out_path = os.path.join(base_dir, "tree_summary.tsv")

    all_rows = []
    all_keys = set()
    tree_col_names = []

    for root, _, _ in os.walk(base_dir):
        nwk_files = glob.glob(os.path.join(root, "seed*.nwk"))
        for nwk_file in nwk_files:
            row = {}
            rel_dir = os.path.relpath(root, PROJECT_ROOT)
            t = dendropy.Tree.get(path=nwk_file, schema="newick")
            tree_params = get_tool_params(nwk_file, config, "tree_sim")
            add_to_ordered_set(tree_col_names, tree_params.keys())
            row.update(tree_params)

            print(f"Processing {rel_dir}...")
            row["tree_path"] = nwk_file
            row["tree_height"] = tree_height(t)
            row["num_leaves"] = num_leaves(t)
            row["colless_index"] = colless_index(t)
            leaf_dist = leaf_distance_to_root(t)
            add_to_ordered_set(tree_col_names, leaf_dist.keys())
            row.update(leaf_dist)
            leaf_steps = leaf_steps_to_root(t)
            add_to_ordered_set(tree_col_names, leaf_steps.keys())
            row.update(leaf_steps)
            pairwise_leaf_dists = pairwise_leaf_step_distances(t)
            add_to_ordered_set(tree_col_names, pairwise_leaf_dists.keys())
            row.update(pairwise_leaf_dists)
            row["max_leaf_dist_to_root"] = max_leaf_distance_to_root(t)
            row["max_leaf_steps_to_root"] = max_leaf_steps_to_root(t)
            
            all_rows.append(row)
            all_keys.update(row.keys())

    column_order = ["tree_path"] + tree_col_names
    remaining_cols = sorted(list(all_keys - set(column_order)))
    column_order += remaining_cols

    write_table(all_rows, column_order, out_path)

if __name__ == "__main__":
    main()
