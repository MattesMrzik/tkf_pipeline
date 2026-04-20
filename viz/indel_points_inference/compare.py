import os
import sys

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import dendropy
from typing import Dict

from viz.indel_points_inference.utils import (
    load_msa,
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

    true_set = set((e.node, e.start, e.end, e.event_type) for e in true_events.events)
    inferred_set = set((e.node, e.start, e.end, e.event_type) for e in inferred_events.events)

    matched_true = true_set & inferred_set
    # unmatched_true = true_set - inferred_set

    nit = sum(1 for e in true_events.events if e.event_type == EventType.INSERTION)
    ndt = sum(1 for e in true_events.events if e.event_type == EventType.DELETION)
    nie = sum(1 for e in inferred_events.events if e.event_type == EventType.INSERTION)
    nde = sum(1 for e in inferred_events.events if e.event_type == EventType.DELETION)

    annotation_agreement = len(matched_true) / len(true_set) if len(true_set) > 0 else 0.0

    if ndt > 0 and nie > 0 and nde > 0:
        true_ratio = nit / ndt
        est_ratio = nie / nde
        indel_ratio = true_ratio / est_ratio
    else:
        indel_ratio = float("nan")

    if nit > 0 or ndt > 0:
        indel_agreement = ((nit - nie) ** 2 + (ndt - nde) ** 2) ** 0.5 / (nit**2 + ndt**2) ** 0.5
    else:
        indel_agreement = 0.0

    return {
        "annotation_agreement": annotation_agreement,
        "indel_ratio": indel_ratio,
        "indel_agreement": indel_agreement,
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
