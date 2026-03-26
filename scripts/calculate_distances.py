import dendropy
import sys
import json
import argparse

def calculate_distances(tree1_path, tree2_path, output_path):
    tree1, tree2 = read_trees(tree1_path, tree2_path)

    rf_distance = dendropy.calculate.treecompare.symmetric_difference(tree1, tree2)
    kf_distance = dendropy.calculate.treecompare.euclidean_distance(tree1, tree2)

    write_results(output_path, rf_distance, kf_distance)

def write_results(output_path, rf_distance, kf_distance):
    results = {
        "robinson_foulds": float(rf_distance),
        "kuhner_felsenstein": float(kf_distance)
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
    parser.add_argument("--tree1", required=True, help="Path to the first tree (Newick)")
    parser.add_argument("--tree2", required=True, help="Path to the second tree (Newick)")
    parser.add_argument("--output", required=True, help="Path to the output JSON file")
    
    args = parser.parse_args()
    calculate_distances(args.tree1, args.tree2, args.output)
