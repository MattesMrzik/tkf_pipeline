import os
import sys
import dendropy
from typing import Dict

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from viz.msa.utils import load_msa
from viz.indel_points_inference.utils import (
    load_tree,
    infer_indels,
    EventType,
)

def compare_indel_annotations(
    tree: dendropy.Tree,
    true_msa: Dict[str, str],
    inferred_msa: Dict[str, str],
) -> Dict[str, float]:
    true_events = infer_indels(true_msa, tree)
    inferred_events = infer_indels(inferred_msa, tree)

    true_set = set(true_events.events)
    inferred_set = set(inferred_events.events)

    matched_true = true_set & inferred_set

    long_nit = true_events.count_by_type(EventType.INSERTION)
    long_ndt = true_events.count_by_type(EventType.DELETION)
    long_nie = inferred_events.count_by_type(EventType.INSERTION)
    long_nde = inferred_events.count_by_type(EventType.DELETION)

    long_annotation_agreement = len(matched_true) / len(true_set) if len(true_set) > 0 else 0.0

    if long_ndt > 0 and long_nie > 0 and long_nde > 0:
        true_ratio = long_nit / long_ndt
        est_ratio = long_nie / long_nde
        long_indel_ratio = true_ratio / est_ratio
    else:
        long_indel_ratio = float("nan")

    if long_nit > 0 or long_ndt > 0:
        nom = ((long_nit - long_nie) ** 2 + (long_ndt - long_nde) ** 2)
        denom = long_nit**2 + long_ndt**2
        long_annotation_agreement = (nom / denom) ** 0.5
    else:
        long_indel_agreement = 0.0

    return {
        "long_annotation_agreement": long_annotation_agreement,
        "long_indel_ratio": long_indel_ratio,
        "long_indel_agreement": long_indel_agreement,
        "long_nit": long_nit,
        "long_ndt": long_ndt,
        "long_nie": long_nie,
        "long_nde": long_nde,
    }


def compare_from_files(
    tree_path: str,
    true_msa_path: str,
    inferred_msa_path: str,
) -> Dict[str, float]:
    tree = load_tree(tree_path)
    true_msa = load_msa(true_msa_path)
    inferred_msa = load_msa(inferred_msa_path)
    return compare_indel_annotations(tree, true_msa, inferred_msa)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Path to tree Newick file")
    parser.add_argument("true_msa", help="Path to true MSA FASTA")
    parser.add_argument("inferred_msa", help="Path to inferred MSA FASTA")
    args = parser.parse_args()

    result = compare_from_files(args.tree, args.true_msa, args.inferred_msa)
    print(result)
