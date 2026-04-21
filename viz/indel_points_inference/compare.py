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
    IndelEvents,
)

def _compute_indel_measures(
    true_events: IndelEvents,
    inferred_events: IndelEvents,
    suffix: str,
) -> Dict[str, float]:
    true_set = set(true_events.events)
    inferred_set = set(inferred_events.events)

    matched_true = true_set & inferred_set

    nit = true_events.count_by_type(EventType.INSERTION)
    ndt = true_events.count_by_type(EventType.DELETION)
    nie = inferred_events.count_by_type(EventType.INSERTION)
    nde = inferred_events.count_by_type(EventType.DELETION)

    # see kimIndelignProbabilisticFramework2007, best 1
    annotation_agreement = len(matched_true) / len(true_set) if len(true_set) > 0 else 0.0

    if ndt > 0 and nie > 0 and nde > 0:
        true_ratio = nit / ndt
        est_ratio = nie / nde
        # see kimIndelignProbabilisticFramework2007, best 1
        indel_ratio = true_ratio / est_ratio
    else:
        indel_ratio = float("nan")

    if nit > 0 or ndt > 0:
        nom = ((nit - nie) ** 2 + (ndt - nde) ** 2)
        denom = nit**2 + ndt**2
        # see kimIndelignProbabilisticFramework2007, best 0
        indel_agreement = (nom / denom) ** 0.5
    else:
        indel_agreement = 0.0

    return {
        f"{suffix}_annotation_agreement": annotation_agreement,
        f"{suffix}_indel_ratio": indel_ratio,
        f"{suffix}_indel_agreement": indel_agreement,
        f"{suffix}_nit": nit,
        f"{suffix}_ndt": ndt,
        f"{suffix}_nie": nie,
        f"{suffix}_nde": nde,
    }


def compare_indel_annotations(
    tree: dendropy.Tree,
    true_msa: Dict[str, str],
    inferred_msa: Dict[str, str],
) -> Dict[str, float]:
    true_events = infer_indels(true_msa, tree)
    inferred_events = infer_indels(inferred_msa, tree)
    long_measures = _compute_indel_measures(true_events, inferred_events, "long")

    true_short = true_events.split_to_single_site()
    inferred_short = inferred_events.split_to_single_site()
    short_measures = _compute_indel_measures(true_short, inferred_short, "short")

    return {**long_measures, **short_measures}


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
