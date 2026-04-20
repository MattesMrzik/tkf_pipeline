import os
import sys


project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from viz.indel_points_inference.compare import compare_from_files, compare_indel_annotations
from viz.indel_points_inference.utils import  load_tree
from viz.msa.utils import load_msa

def test_compare_indel_annotations():
    base_dir = os.path.join(os.path.dirname(__file__), "data")
    tree_path = os.path.join(base_dir, "tree.nwk")
    true_msa_path = os.path.join(base_dir, "true_msa.fasta")
    inferred_msa_path = os.path.join(base_dir, "inferred_msa.fasta")

    result = compare_from_files(tree_path, true_msa_path, inferred_msa_path)

    print(f"Annotation agreement: {result['long_annotation_agreement']}")
    print(f"Indel ratio: {result['long_indel_ratio']}")
    print(f"Indel agreement: {result['long_indel_agreement']}")
    print(f"long_nit={result['long_nit']}, long_ndt={result['long_ndt']}, long_nie={result['long_nie']}, long_nde={result['long_nde']}")

    assert result["long_annotation_agreement"] == 0.0
    assert result["long_indel_ratio"] == 0.75
    assert result["long_indel_agreement"] > 0.3
    assert result["long_nit"] == 3
    assert result["long_ndt"] == 1
    assert result["long_nie"] == 4
    assert result["long_nde"] == 1

    tree = load_tree(tree_path)
    true_msa = load_msa(true_msa_path)
    inferred_msa = load_msa(inferred_msa_path)

    result2 = compare_indel_annotations(tree, true_msa, inferred_msa)
    assert result == result2

    print("All tests passed!")

if __name__ == "__main__":
    test_compare_indel_annotations()
