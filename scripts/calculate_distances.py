import dendropy
import sys
import json
import argparse

def calculate_distances(true_tree_path, final_tree_path, start_tree_path, output_path):
    true_tree, final_tree = read_trees(true_tree_path, final_tree_path)
    
    # Calculate distances for final tree
    final_rf = dendropy.calculate.treecompare.symmetric_difference(true_tree, final_tree)
    final_kf = dendropy.calculate.treecompare.euclidean_distance(true_tree, final_tree)

    # Calculate distances for start tree
    # We need to re-read or re-reconcile because TaxonNamespace needs to be consistent
    # Let's just read start_tree and reconcile it with the same TaxonNamespace
    start_tree_obj = dendropy.Tree.get(
        path=start_tree_path,
        schema="newick",
        preserve_underscores=True,
        suppress_internal_node_taxa=True,
        suppress_leaf_node_taxa=True,
        taxon_namespace=true_tree.taxon_namespace
    )
    
    def reconcile_taxa(tree, namespace):
        for leaf in tree.leaf_node_iter():
            if leaf.label:
                taxon = namespace.require_taxon(label=leaf.label)
                leaf.taxon = taxon
    
    reconcile_taxa(start_tree_obj, true_tree.taxon_namespace)
    start_tree_obj.is_rooted = False
    
    start_rf = dendropy.calculate.treecompare.symmetric_difference(true_tree, start_tree_obj)
    start_kf = dendropy.calculate.treecompare.euclidean_distance(true_tree, start_tree_obj)

    results = {
        "robinson_foulds": float(final_rf),
        "kuhner_felsenstein": float(final_kf),
        "start_robinson_foulds": float(start_rf),
        "start_kuhner_felsenstein": float(start_kf)
    }
    
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=4)

def read_trees(tree1_path, tree2_path):
    try:
        # We use suppress_leaf_node_taxa to avoid "duplicate taxon labels" errors,
        # then we map them to a unified TaxonNamespace for comparison.
        tree1 = dendropy.Tree.get(
            path=tree1_path,
            schema="newick",
            preserve_underscores=True,
            suppress_internal_node_taxa=True,
            suppress_leaf_node_taxa=True
        )
        tree2 = dendropy.Tree.get(
            path=tree2_path,
            schema="newick",
            preserve_underscores=True,
            suppress_internal_node_taxa=True,
            suppress_leaf_node_taxa=True
        )
    except Exception as e:
        print(f"Error loading trees: {e}", file=sys.stderr)
        sys.exit(1)

    # Re-associate labels with a unified TaxonNamespace
    tns = dendropy.TaxonNamespace()
    
    # Dendropy does not automatically map labels to Taxon objects if they were suppressed.
    # We manually populate the namespace and link the nodes.
    def reconcile_taxa(tree, namespace):
        for leaf in tree.leaf_node_iter():
            if leaf.label:
                taxon = namespace.require_taxon(label=leaf.label)
                leaf.taxon = taxon

    reconcile_taxa(tree1, tns)
    reconcile_taxa(tree2, tns)
    
    # Set the namespace for both trees
    tree1.taxon_namespace = tns
    tree2.taxon_namespace = tns

    tree1.is_rooted = False
    tree2.is_rooted = False

    return (tree1, tree2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RF and KF distances between two trees.")
    parser.add_argument("--true-tree", required=True, help="Path to the true tree (Newick)")
    parser.add_argument("--final-tree", required=True, help="Path to the final inferred tree (Newick)")
    parser.add_argument("--start-tree", required=True, help="Path to the starting tree (NJ) (Newick)")
    parser.add_argument("--output", required=True, help="Path to the output JSON file")
    
    args = parser.parse_args()
    calculate_distances(args.true_tree, args.final_tree, args.start_tree, args.output)
