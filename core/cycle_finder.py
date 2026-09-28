"""
core/cycle_finder.py
====================
Automatic detection and analysis of cyclic transitions in an isotope network.

Algorithm:
  1. Builds a directed graph: nodes are isotopes, edges are decays and reactions.
  2. Tarjan's SCC finds strongly connected components (groups of mutually
     reachable isotopes). Each SCC with >= 2 vertices contains >= 1 cycle.
  3. DFS within each SCC enumerates all elementary cycles.
  4. Canonicalization (rotate to lexicographically minimum vertex) removes duplicates.

Input: list of Isotope objects after IsotopeBuilder.link_products().

Usage example:
    from core import CycleAnalyzer
    analyzer = CycleAnalyzer(isotopes)     # after build_range / link_products
    cycles   = analyzer.find_all_cycles()
    analyzer.print_report()
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

from .entities import Isotope


# ─────────────────────────────── Data structures ──────────────────────────

@dataclass
class CycleStep:
    """
    One transition in the cycle: from 'isotope' via 'label' to the next isotope.
    label — human-readable reaction or decay label.
    kind  — "decay" | "reaction"
    mt    — ENDF MT number for reactions; None for decays.
    """
    isotope: Isotope
    label:   str
    kind:    str
    mt:      Optional[int] = None

    def __str__(self) -> str:
        return f"{self.isotope.name} --[{self.label}]-->"


@dataclass
class IsotopeCycle:
    """
    Elementary cycle: steps[0], steps[1], ..., steps[n-1],
    where steps[n-1] leads back to steps[0].isotope.
    """
    steps: List[CycleStep]

    # ── properties ──

    @property
    def length(self) -> int:
        return len(self.steps)

    @property
    def isotope_names(self) -> List[str]:
        return [s.isotope.name for s in self.steps]

    @property
    def is_pure_decay(self) -> bool:
        return all(s.kind == "decay" for s in self.steps)

    @property
    def is_pure_reaction(self) -> bool:
        return all(s.kind == "reaction" for s in self.steps)

    @property
    def is_mixed(self) -> bool:
        kinds = {s.kind for s in self.steps}
        return len(kinds) > 1

    @property
    def n_decays(self) -> int:
        return sum(1 for s in self.steps if s.kind == "decay")

    @property
    def n_reactions(self) -> int:
        return sum(1 for s in self.steps if s.kind == "reaction")

    def __len__(self) -> int:
        return self.length

    def __str__(self) -> str:
        parts = [str(s) for s in self.steps]
        return " ".join(parts) + f" {self.steps[0].isotope.name}"


@dataclass
class SCC:
    """
    Strongly Connected Component — group of mutually reachable isotopes.
    The cycles field is populated after find_all_cycles() is called.
    """
    isotopes: List[Isotope]
    cycles:   List[IsotopeCycle] = field(default_factory=list)

    @property
    def names(self) -> List[str]:
        return sorted(iso.name for iso in self.isotopes)

    def __len__(self) -> int:
        return len(self.isotopes)

    def __repr__(self) -> str:
        return f"SCC({self.names})"


# ─────────────────────────────── Main class ────────────────────────────────

class CycleAnalyzer:
    """
    Analyzer for cyclic transitions in an isotope network.

    Parameters
    ----------
    isotopes : list of Isotope objects with product references set
               (after IsotopeBuilder.link_products or build_range).

    Example
    -------
    analyzer = CycleAnalyzer(isotopes)
    cycles   = analyzer.find_all_cycles()
    analyzer.print_report()
    """

    _MT_NAMES: Dict[int, str] = {
        102: "(n,g)", 103: "(n,p)", 107: "(n,a)",
        16:  "(n,2n)", 17: "(n,3n)", 18: "fission",
    }

    def __init__(self, isotopes: List[Isotope]):
        self._isotopes: List[Isotope] = list(isotopes)
        self._iso_set:  Set[Isotope]  = set(isotopes)
        # adjacency graph: iso -> [(target, label, kind, mt), ...]
        self._adj: Dict[Isotope, List[Tuple[Isotope, str, str, Optional[int]]]] = {}
        self._build_adj()

    # ── graph construction ──

    def _build_adj(self) -> None:
        for iso in self._isotopes:
            edges: List[Tuple[Isotope, str, str, Optional[int]]] = []

            for d in iso.getListOfDecays():
                if d.product is not None and d.product in self._iso_set:
                    lbl = d.mode
                    if d.branch < 0.999:
                        lbl += f" [{d.branch * 100:.2g}%]"
                    edges.append((d.product, lbl, "decay", None))

            for r in iso.getListOfReactions():
                if r.product is not None and r.product in self._iso_set:
                    mt_name = self._MT_NAMES.get(r.mt, f"MT{r.mt}")
                    lbl = f"{mt_name} s={r.sigma_barn:.3g}b"
                    edges.append((r.product, lbl, "reaction", r.mt))

            self._adj[iso] = edges

    def _count_edges(self) -> int:
        return sum(len(v) for v in self._adj.values())

    # ── Tarjan's SCC (recursive version, no depth restrictions) ──

    def find_sccs(self) -> List[SCC]:
        """
        Tarjan's algorithm — O(V+E).
        Returns only non-trivial SCCs (>= 2 vertices or a self-loop).
        """
        idx:      Dict[Isotope, int] = {}
        low:      Dict[Isotope, int] = {}
        on_stack: Set[Isotope]       = set()
        stack:    List[Isotope]      = []
        counter = [0]
        result:   List[SCC]          = []

        def visit(v: Isotope) -> None:
            idx[v] = low[v] = counter[0]
            counter[0] += 1
            stack.append(v)
            on_stack.add(v)

            for (w, _, _, _) in self._adj.get(v, []):
                if w not in idx:
                    visit(w)
                    low[v] = min(low[v], low[w])
                elif w in on_stack:
                    low[v] = min(low[v], idx[w])

            if low[v] == idx[v]:
                comp: List[Isotope] = []
                while True:
                    w = stack.pop()
                    on_stack.discard(w)
                    comp.append(w)
                    if w is v:
                        break
                if len(comp) > 1:
                    result.append(SCC(comp))

        for v in self._isotopes:
            if v not in idx:
                visit(v)

        return result

    # ── enumerate all elementary cycles ──

    def find_all_cycles(self) -> List[IsotopeCycle]:
        """
        Enumerates all elementary cycles in the network.

        Method:
          - For each non-trivial SCC, runs DFS from each vertex.
          - Canonical representation (rotate to min-name) removes duplicates.

        Returns a list of IsotopeCycle sorted by length.
        """
        sccs = self.find_sccs()
        all_cycles: List[IsotopeCycle] = []
        seen_keys:  Set[tuple]         = set()

        for scc in sccs:
            scc_set = set(scc.isotopes)
            # subgraph — only edges within the SCC
            sub: Dict[Isotope, List] = {
                v: [(w, lbl, kind, mt)
                    for (w, lbl, kind, mt) in self._adj.get(v, [])
                    if w in scc_set]
                for v in scc.isotopes
            }

            for start in scc.isotopes:
                for cycle in self._dfs_from(start, scc_set, sub):
                    key = self._canon_key(cycle)
                    if key not in seen_keys:
                        seen_keys.add(key)
                        all_cycles.append(cycle)
                        scc.cycles.append(cycle)

        all_cycles.sort(key=lambda c: c.length)
        return all_cycles

    def _dfs_from(
        self,
        start:   Isotope,
        scc_set: Set[Isotope],
        sub:     Dict[Isotope, List],
    ) -> List[IsotopeCycle]:
        """DFS from vertex start — finds all elementary cycles through it."""
        results:    List[IsotopeCycle] = []
        path_nodes: Set[Isotope]       = {start}

        def dfs(v: Isotope, steps: List[CycleStep]) -> None:
            for (w, lbl, kind, mt) in sub.get(v, []):
                step = CycleStep(v, lbl, kind, mt)
                if w is start and steps:
                    # Cycle closed: steps describe path from start to v,
                    # step is the transition v -> start
                    results.append(IsotopeCycle(steps + [step]))
                elif w not in path_nodes:
                    path_nodes.add(w)
                    dfs(w, steps + [step])
                    path_nodes.discard(w)

        dfs(start, [])
        return results

    @staticmethod
    def _canon_key(cycle: IsotopeCycle) -> tuple:
        """
        Canonical cycle representation: rotates steps so that the isotope
        with the lexicographically smallest name comes first.
        This deduplicates cycles found from different starting vertices.
        """
        pairs = [(s.isotope.name, s.label) for s in cycle.steps]
        mi = min(range(len(pairs)), key=lambda i: pairs[i][0])
        rotated = pairs[mi:] + pairs[:mi]
        return tuple(rotated)

    # ── queries without full enumeration ──

    def has_cycles(self) -> bool:
        """Quick check for at least one cycle (Tarjan only)."""
        return bool(self.find_sccs())

    def isotopes_in_cycles(self) -> List[Isotope]:
        """List of isotopes participating in at least one cycle."""
        result: List[Isotope] = []
        for scc in self.find_sccs():
            result.extend(scc.isotopes)
        return result

    def cycles_through(self, isotope: Isotope) -> List[IsotopeCycle]:
        """All cycles passing through the given isotope."""
        return [c for c in self.find_all_cycles()
                if any(s.isotope is isotope for s in c.steps)]

    # ── report ──

    def print_report(
        self,
        cycles:   Optional[List[IsotopeCycle]] = None,
        max_show: int = 30,
    ) -> None:
        """
        Prints a full report: SCC groups, cycle summary, details.

        Parameters
        ----------
        cycles   : if None — computed automatically
        max_show : maximum number of cycles shown per category
        """
        if cycles is None:
            cycles = self.find_all_cycles()

        sccs = self.find_sccs()
        sep = "=" * 72

        print(f"\n{sep}")
        print(f"  CYCLE ANALYSIS  "
              f"({len(self._isotopes)} isotopes, {self._count_edges()} edges)")
        print(sep)

        # ── SCC ──
        print(f"\n  Strongly Connected Components ({len(sccs)} non-trivial):")
        if not sccs:
            print("    None — network is acyclic (DAG).")
        for scc in sccs:
            print(f"    SCC [{len(scc)} isotopes]: {', '.join(scc.names)}")

        # ── summary ──
        print(f"\n  Elementary cycles found: {len(cycles)}")
        if not cycles:
            print("    None — no reaction loops detected.")
            print(f"\n{sep}")
            return

        pure_dec = [c for c in cycles if c.is_pure_decay]
        pure_rxn = [c for c in cycles if c.is_pure_reaction]
        mixed    = [c for c in cycles if c.is_mixed]

        def _show_group(title: str, group: List[IsotopeCycle]) -> None:
            if not group:
                return
            print(f"\n  -- {title} ({len(group)}) --")
            for i, c in enumerate(group[:max_show], 1):
                type_tag = ""
                if c.is_mixed:
                    type_tag = f"  [{c.n_decays}dec + {c.n_reactions}rxn]"
                print(f"    [{i:2d}] L={c.length}{type_tag}")
                print(f"          {c}")
            if len(group) > max_show:
                print(f"          ... +{len(group) - max_show} more")

        _show_group("Mixed cycles (decay + reaction)", mixed)
        _show_group("Pure-decay cycles", pure_dec)
        _show_group("Pure-reaction cycles", pure_rxn)

        # ── participating isotopes ──
        print(f"\n  Isotopes participating in cycles:")
        in_cycle = {s.isotope.name for c in cycles for s in c.steps}
        not_in   = [i.name for i in self._isotopes if i.name not in in_cycle]
        print(f"    In cycles : {', '.join(sorted(in_cycle))}")
        if not_in:
            print(f"    Not in any: {', '.join(not_in)}")

        print(f"\n{sep}\n")
