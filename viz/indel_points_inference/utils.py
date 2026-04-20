import dendropy
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Set, Optional
from intervaltree import IntervalTree


class EventType(Enum):
    INSERTION = "insertion"
    DELETION = "deletion"

# This annotation is immutable to ensure that events can be safely used in sets and as dictionary keys.
@dataclass(frozen=True)
class IndelEvent:
    node: str
    start: int
    end: int
    event_type: EventType

    def overlaps_column(self, col: int) -> bool:
        return self.start <= col < self.end


class IndelEvents:
    def __init__(self, events: Optional[List[IndelEvent]] = None):
        self.events: List[IndelEvent] = events or []
        self._by_node: Dict[str, List[IndelEvent]] = {}
        self._interval_tree_for_events = IntervalTree()
        self._build_indices()

    def _build_indices(self):
        self._by_node.clear()
        self._interval_tree_for_events = IntervalTree()
        for event in self.events:
            if event.node not in self._by_node:
                self._by_node[event.node] = []
            self._by_node[event.node].append(event)
            self._interval_tree_for_events.addi(event.start, event.end, event)

    def add(self, event: IndelEvent):
        self.events.append(event)
        if event.node not in self._by_node:
            self._by_node[event.node] = []
        self._by_node[event.node].append(event)
        self._interval_tree_for_events.addi(event.start, event.end, event)

    def get_by_node(self, node: str) -> List[IndelEvent]:
        return self._by_node.get(node, [])

    def get_by_column(self, col: int) -> List[IndelEvent]:
        intervals = self._interval_tree_for_events.overlap(col, col + 1)
        return [interval.data for interval in intervals]

    def get_columns_in_region(self, start: int, end: int) -> Set[int]:
        intervals = self._interval_tree_for_events.overlap(start, end)
        columns = set()
        for interval in intervals:
            ev = interval.data
            columns.update(range(max(ev.start, start), min(ev.end, end)))
        return columns

    def count_by_type(self, event_type: EventType) -> int:
        return sum(1 for e in self.events if e.event_type == event_type)


def infer_indels(msa: Dict[str, str], tree: dendropy.Tree) -> IndelEvents:
    events = IndelEvents()
    node_to_seq = {node.label: msa[node.label] for node in tree.preorder_node_iter() if node.label}

    msa_len = len(next(iter(msa.values())))
    for node in tree.preorder_node_iter():
        if node.label is None:
            continue
        if node.parent_node is None:
            continue

        parent_label = node.parent_node.label
        if parent_label is None or parent_label not in node_to_seq:
            continue

        child_seq = node_to_seq[node.label]
        parent_seq = node_to_seq[parent_label]

        i = 0
        while i < msa_len:
            if parent_seq[i] == "-" and child_seq[i] != "-":
                start = i
                while i < msa_len and parent_seq[i] == "-" and child_seq[i] != "-":
                    i += 1
                end = i
                events.add(IndelEvent(node=node.label, start=start, end=end, event_type=EventType.INSERTION))
            elif parent_seq[i] != "-" and child_seq[i] == "-":
                start = i
                while i < msa_len and parent_seq[i] != "-" and child_seq[i] == "-":
                    i += 1
                end = i
                events.add(IndelEvent(node=node.label, start=start, end=end, event_type=EventType.DELETION))
            else:
                i += 1

    return events



def load_tree(newick_path: str) -> dendropy.Tree:
    return dendropy.Tree.get(
        path=newick_path,
        schema="newick",
        preserve_underscores=True,
        suppress_internal_node_taxa=True,
        suppress_leaf_node_taxa=True,
    )
