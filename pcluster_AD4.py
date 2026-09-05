#!/usr/bin/env python3
"""
dockcluster3d_dlg.py

A production-quality script for AutoDock4 DLG cluster analysis.

Reads an AutoDock4 .dlg docking log file, extracts ligand coordinates for each
docking conformation, calculates mass-weighted center-of-mass (COM), clusters
the conformations by COM-to-COM distance (transitive, cutoff 1.0 Å by default),
extracts six energy parameters per conformation, and outputs a cluster.csv
(one row per cluster, representative by lowest binding energy) and a detailed
cluster.txt log.

Usage examples:
    python dockcluster3d_dlg.py --input docking.dlg
    python dockcluster3d_dlg.py --input docking.dlg --output ./clustered --cutoff 0.8

    # Show the script version
    python dockcluster3d_dlg.py --version

Dependencies: None (pure Python 3.6+ standard library).

Version: 1.0

Changelog:
    1.0 - Initial refinement pass. Removed dead/no-op code (an abandoned
          first energy-parsing loop and a redundant duplicate branch for the
          Inhibition Constant that the generic parser already handled),
          removed an unused NamedTuple and an unused import, fixed a
          formatting bug where very small/large-magnitude values (e.g.
          Estimated Inhibition Constant, which is typically 1e-6 to 1e-9 M)
          were rounded to "0.0000" by a fixed 4-decimal format in
          cluster.txt, fixed a bug where missing energy values printed the
          literal text "None" in cluster.csv instead of "N/A", wired the
          parser's returned stats into the console summary instead of a
          hardcoded "0", and added --version / __version__.
"""

import os
import sys
import re
import math
import argparse
import warnings
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, NamedTuple

__version__ = "1.0"

# ----------------------------------------------------------------------
# Atomic masses (in g/mol, but any consistent unit works as they cancel)
ATOMIC_MASS: Dict[str, float] = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "P": 30.974,
    "S": 32.065,
    "Cl": 35.45,
    "Br": 79.904,
    "I": 126.904,
    "B": 10.81,
    "Si": 28.085,
    # Additional common elements:
    "Na": 22.990,
    "Mg": 24.305,
    "K": 39.098,
    "Ca": 40.078,
    "Fe": 55.845,
    "Zn": 65.38,
    "Se": 78.971,
    "Mn": 54.938,
    "Cu": 63.546,
    "Co": 58.933,
}

ENERGY_KEYS = [
    'Estimated Free Energy of Binding',
    'Estimated Inhibition Constant',
    'Final Intermolecular Energy',
    'Final Total Internal Energy',
    'Torsional Free Energy',
    "Unbound System's Energy",
]


class Atom(NamedTuple):
    """Representation of a single atom from a DOCKED: record."""
    x: float
    y: float
    z: float
    element: str
    mass: float


class Conformation(NamedTuple):
    """All data for one docking conformation."""
    conf_id: int                 # internal consecutive ID (1-based)
    atoms: List[Atom]
    com: Tuple[float, float, float]   # (x, y, z) mass-weighted
    energies: Dict[str, Optional[float]]  # six energy keys


# ----------------------------------------------------------------------
# DLG parsing
def parse_pdb_coordinate_line(line: str) -> Optional[Atom]:
    """
    Parse a line that looks like a PDB/PDBQT ATOM/HETATM record (including a
    DOCKED: prefix). Coordinates are read from fixed columns 31-38, 39-46,
    47-54 (1-based), matching the standard PDB/AutoDock layout.

    Returns an Atom, or None if the line cannot be parsed.
    """
    if not line.startswith(('ATOM', 'HETATM')):
        # The line may carry a "DOCKED:" prefix before the ATOM/HETATM record.
        if 'DOCKED:' in line:
            idx = line.find('ATOM')
            if idx == -1:
                idx = line.find('HETATM')
            if idx != -1:
                line = line[idx:]
            else:
                return None
        else:
            return None

    if len(line) < 54:
        return None

    try:
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
    except ValueError:
        return None

    # Determine element symbol.
    element = None
    # 1) Try the element column (73-76), if present.
    if len(line) >= 76:
        elem_col = line[72:76].strip()
        if elem_col:
            element = elem_col.rstrip('0123456789')
            if element not in ATOMIC_MASS:
                element = None
    # 2) Fall back to inferring from the atom name (columns 13-16).
    if element is None and len(line) >= 16:
        atom_name = line[12:16].strip()
        letters = ''.join(ch for ch in atom_name if ch.isalpha())
        if letters:
            if letters[0] in ATOMIC_MASS:
                element = letters[0]
            elif len(letters) >= 2 and letters[:2] in ATOMIC_MASS:
                element = letters[:2]
            elif len(letters) >= 2 and letters[1] in ATOMIC_MASS:
                element = letters[1]
            else:
                candidate = letters.upper()
                if candidate in ATOMIC_MASS:
                    element = candidate
                else:
                    candidate = letters.capitalize()
                    if candidate in ATOMIC_MASS:
                        element = candidate
    # 3) Last resort: the residue name (columns 18-20) sometimes is the element.
    if element is None and len(line) >= 20:
        res_name = line[17:20].strip()
        if res_name in ATOMIC_MASS:
            element = res_name
    if element is None:
        return None

    mass = ATOMIC_MASS.get(element)
    if mass is None:
        return None
    return Atom(x=x, y=y, z=z, element=element, mass=mass)


def _extract_energy(line: str, key: str) -> Optional[float]:
    """Find `key ... = <number>` in a line (case-insensitive) and return the number."""
    pattern = re.compile(re.escape(key) + r'.*?=\s*([\d.eE+-]+)', re.IGNORECASE)
    match = pattern.search(line)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None


def extract_conformations_from_dlg(filepath: str) -> Tuple[List[Conformation], Dict]:
    """
    Parse the DLG file and extract all docking conformations with their
    coordinates and energies.

    Returns:
        - list of Conformation objects (ordered as they appear)
        - dict with statistics: total_parsed, skipped, warnings
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    stats_warnings: List[str] = []

    # ------------------------------------------------------------------
    # Pass 1: find DOCKED: blocks and extract atoms. A block starts at the
    # first "DOCKED: ... ATOM/HETATM" line and ends at a "DOCKED: ... END"
    # line (matches END, ENDMDL, ENDROOT, etc.).
    conformations = []
    current_atoms: List[Atom] = []
    in_conformation = False
    conf_counter = 0

    for line in lines:
        if 'DOCKED:' in line and ('ATOM' in line or 'HETATM' in line):
            if not in_conformation:
                in_conformation = True
                current_atoms = []
            atom = parse_pdb_coordinate_line(line)
            if atom is not None:
                current_atoms.append(atom)
        elif 'DOCKED:' in line and 'END' in line:
            if in_conformation and current_atoms:
                conf_counter += 1
                conformations.append({
                    'conf_id': conf_counter,
                    'atoms': current_atoms,
                    'com': compute_com(current_atoms),
                    'energies': {},
                })
            in_conformation = False
            current_atoms = []

    # Close a dangling block if the file didn't end with an explicit END line.
    if in_conformation and current_atoms:
        conf_counter += 1
        conformations.append({
            'conf_id': conf_counter,
            'atoms': current_atoms,
            'com': compute_com(current_atoms),
            'energies': {},
        })

    # ------------------------------------------------------------------
    # Pass 2: parse the per-run energy summary lines, e.g.:
    #   Run 1: Estimated Free Energy of Binding = -9.21 kcal/mol
    #   Run 1: Estimated Inhibition Constant, Ki = 2.45e-06 M
    #   Run 1: Final Intermolecular Energy = -10.32 kcal/mol
    #   ...
    run_energies: Dict[int, Dict[str, float]] = defaultdict(dict)
    run_num = None
    for line in lines:
        run_match = re.search(r'Run\s*(\d+):', line)
        if run_match:
            run_num = int(run_match.group(1))
        if run_num is not None:
            for key in ENERGY_KEYS:
                val = _extract_energy(line, key)
                if val is not None:
                    run_energies[run_num][key] = val

    # ------------------------------------------------------------------
    # Associate energies with conformations by matching run number to
    # conformation order (conformation i <-> run i+1).
    num_confs = len(conformations)
    num_runs = len(run_energies)
    if num_confs != num_runs:
        msg = (f"Number of conformations ({num_confs}) does not match number "
               f"of energy runs ({num_runs}). Will try to align by order; "
               f"may be incomplete.")
        warnings.warn(msg)
        stats_warnings.append(msg)

    for i, conf_data in enumerate(conformations):
        run_num = i + 1
        if run_num in run_energies:
            conf_data['energies'] = run_energies[run_num]
        else:
            msg = f"Conformation {conf_data['conf_id']} has no energy data."
            warnings.warn(msg)
            stats_warnings.append(msg)
            conf_data['energies'] = {}

    # Build the final Conformation objects, filling any missing energy key
    # with None so downstream code can rely on all six keys being present.
    result = []
    for data in conformations:
        energies = data['energies']
        for key in ENERGY_KEYS:
            energies.setdefault(key, None)
        result.append(Conformation(
            conf_id=data['conf_id'],
            atoms=data['atoms'],
            com=data['com'],
            energies=energies,
        ))

    stats = {
        'total_parsed': len(result),
        'skipped': 0,
        'warnings': stats_warnings,
    }
    return result, stats


def compute_com(atoms: List[Atom]) -> Tuple[float, float, float]:
    """Compute mass-weighted center of mass from an atom list."""
    total_mass = 0.0
    cx = cy = cz = 0.0
    for atom in atoms:
        m = atom.mass
        total_mass += m
        cx += m * atom.x
        cy += m * atom.y
        cz += m * atom.z
    if total_mass == 0.0:
        return (0.0, 0.0, 0.0)
    return (cx / total_mass, cy / total_mass, cz / total_mass)


# ----------------------------------------------------------------------
# Clustering (Union-Find)
class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x: int, y: int) -> None:
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return
        if self.rank[rx] < self.rank[ry]:
            self.parent[rx] = ry
        elif self.rank[rx] > self.rank[ry]:
            self.parent[ry] = rx
        else:
            self.parent[ry] = rx
            self.rank[rx] += 1


def cluster_conformations(conformations: List[Conformation], cutoff: float) -> List[List[int]]:
    """
    Cluster conformations by COM distances using transitive closure.
    Returns a list of clusters, each a list of indices into `conformations`.
    """
    n = len(conformations)
    uf = UnionFind(n)
    for i in range(n):
        com_i = conformations[i].com
        for j in range(i + 1, n):
            com_j = conformations[j].com
            dist = math.sqrt((com_i[0] - com_j[0]) ** 2 +
                              (com_i[1] - com_j[1]) ** 2 +
                              (com_i[2] - com_j[2]) ** 2)
            if dist <= cutoff:
                uf.union(i, j)
    clusters_dict: Dict[int, List[int]] = defaultdict(list)
    for idx in range(n):
        root = uf.find(idx)
        clusters_dict[root].append(idx)
    return list(clusters_dict.values())


def cluster_com(conformations: List[Conformation], indices: List[int]) -> Tuple[float, float, float]:
    """Compute the mean of the ligand COMs in the cluster."""
    if not indices:
        return (0.0, 0.0, 0.0)
    sum_x = sum_y = sum_z = 0.0
    for idx in indices:
        x, y, z = conformations[idx].com
        sum_x += x
        sum_y += y
        sum_z += z
    n = len(indices)
    return (sum_x / n, sum_y / n, sum_z / n)


def order_clusters(clusters: List[List[int]], conformations: List[Conformation]) -> List[List[int]]:
    """Sort clusters by cluster COM (X, then Y, then Z)."""
    coms = [cluster_com(conformations, cluster) for cluster in clusters]
    pairs = sorted(zip(coms, clusters), key=lambda p: (p[0][0], p[0][1], p[0][2]))
    return [p[1] for p in pairs]


def select_representative(conf_ids: List[int], conformations: List[Conformation]) -> int:
    """
    Choose the representative conformation as the one with the lowest
    Estimated Free Energy of Binding. Returns the conformation ID.
    """
    by_id = {c.conf_id: c for c in conformations}
    best_conf_id = None
    best_energy = float('inf')
    for cid in conf_ids:
        conf = by_id[cid]
        energy = conf.energies.get('Estimated Free Energy of Binding')
        if energy is not None and energy < best_energy:
            best_energy = energy
            best_conf_id = cid
    if best_conf_id is None:
        best_conf_id = conf_ids[0]
    return best_conf_id


def _format_energy(key: str, val: Optional[float]) -> str:
    """
    Format an energy value for the text log. Estimated Inhibition Constant
    (Ki) values are typically 1e-6 to 1e-9 M, so a fixed 4-decimal format
    would print "0.0000" for essentially all of them; use scientific
    notation for that field and fixed-point for the rest.
    """
    if val is None:
        return "N/A"
    if key == 'Estimated Inhibition Constant':
        return f"{val:.3e}"
    return f"{val:.4f}"


# ----------------------------------------------------------------------
# Output writing
def write_cluster_txt(conformations: List[Conformation],
                       clusters: List[List[int]],
                       output_dir: str,
                       cutoff: float,
                       input_file: str) -> None:
    """Write the detailed cluster.txt log."""
    path = os.path.join(output_dir, 'cluster.txt')
    with open(path, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("DockCluster3D DLG Analysis - Cluster Log\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Script version: {__version__}\n")
        f.write(f"Input DLG: {input_file}\n")
        f.write(f"Cutoff: {cutoff:.3f} \u00c5\n")
        f.write(f"Total conformations: {len(conformations)}\n")
        f.write(f"Number of clusters: {len(clusters)}\n\n")

        for i, cluster_indices in enumerate(clusters, start=1):
            conf_ids = sorted(conformations[idx].conf_id for idx in cluster_indices)
            com = cluster_com(conformations, cluster_indices)
            rep_id = select_representative(conf_ids, conformations)
            rep_conf = next(c for c in conformations if c.conf_id == rep_id)

            f.write("=" * 60 + "\n")
            f.write(f"CLUSTER {i}\n")
            f.write("=" * 60 + "\n")
            f.write(f"Number of conformations: {len(conf_ids)}\n")
            f.write(f"Cluster COM: X = {com[0]:.3f}, Y = {com[1]:.3f}, Z = {com[2]:.3f}\n")
            f.write(f"Representative conformation: {rep_id}\n")
            f.write("Representative energies:\n")
            for key in ENERGY_KEYS:
                f.write(f"  {key}: {_format_energy(key, rep_conf.energies.get(key))}\n")
            f.write("\nConformation IDs in cluster:\n")
            for cid in conf_ids:
                conf = next(c for c in conformations if c.conf_id == cid)
                com_lig = conf.com
                f.write(f"  {cid:3d}  COM: ({com_lig[0]:.3f}, {com_lig[1]:.3f}, {com_lig[2]:.3f})\n")
            f.write("\n")


def write_cluster_csv(conformations: List[Conformation],
                       clusters: List[List[int]],
                       output_dir: str) -> None:
    """Write cluster.csv with one row per cluster."""
    path = os.path.join(output_dir, 'cluster.csv')
    with open(path, 'w') as f:
        header = ["Cluster", "Conformations"] + ENERGY_KEYS
        f.write(",".join(header) + "\n")

        for i, cluster_indices in enumerate(clusters, start=1):
            conf_ids = sorted(conformations[idx].conf_id for idx in cluster_indices)
            rep_id = select_representative(conf_ids, conformations)
            rep_conf = next(c for c in conformations if c.conf_id == rep_id)

            row = [str(i), '"' + ",".join(str(cid) for cid in conf_ids) + '"']
            for key in ENERGY_KEYS:
                val = rep_conf.energies.get(key)
                row.append("N/A" if val is None else str(val))
            f.write(",".join(row) + "\n")


# ----------------------------------------------------------------------
# Main CLI
def main():
    parser = argparse.ArgumentParser(
        description="Cluster AutoDock4 DLG conformations by center-of-mass and generate cluster.csv.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--version", action="version",
                         version=f"%(prog)s {__version__}",
                         help="Show the script version and exit")
    parser.add_argument("--input", required=True,
                         help="AutoDock4 .dlg file to parse")
    parser.add_argument("--output", default=None,
                         help="Output directory (default: same as input file location)")
    parser.add_argument("--cutoff", type=float, default=1.0,
                         help="Clustering cutoff in \u00c5 (default: 1.0)")

    args = parser.parse_args()

    input_file = args.input
    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        sys.exit(1)

    output_dir = args.output if args.output is not None else (os.path.dirname(input_file) or ".")
    os.makedirs(output_dir, exist_ok=True)

    cutoff = args.cutoff
    if cutoff < 0:
        print("Error: Cutoff must be non-negative.", file=sys.stderr)
        sys.exit(1)

    # Parse DLG
    try:
        conformations, stats = extract_conformations_from_dlg(input_file)
    except Exception as e:
        print(f"Error parsing DLG: {e}", file=sys.stderr)
        sys.exit(1)

    if not conformations:
        print("No docking conformations found in DLG.", file=sys.stderr)
        with open(os.path.join(output_dir, 'cluster.txt'), 'w') as f:
            f.write("No conformations found.\n")
        with open(os.path.join(output_dir, 'cluster.csv'), 'w') as f:
            f.write(",".join(["Cluster", "Conformations"] + ENERGY_KEYS) + "\n")
        sys.exit(0)

    # Cluster
    raw_clusters = cluster_conformations(conformations, cutoff)
    ordered_clusters = order_clusters(raw_clusters, conformations)

    # Write outputs
    write_cluster_txt(conformations, ordered_clusters, output_dir, cutoff, input_file)
    write_cluster_csv(conformations, ordered_clusters, output_dir)

    # Console summary
    print("=" * 60)
    print("DockCluster3D DLG Analysis Complete")
    print("=" * 60)
    print(f"Script version: {__version__}")
    print(f"Input DLG: {input_file}")
    print(f"Conformations parsed: {stats['total_parsed']}")
    print(f"Conformations skipped: {stats['skipped']}")
    print(f"Warnings encountered: {len(stats['warnings'])}")
    print(f"Clusters generated: {len(ordered_clusters)}")
    print(f"Clustering cutoff: {cutoff:.3f} \u00c5")
    print(f"Output directory: {output_dir}/")
    print(f"CSV: {os.path.join(output_dir, 'cluster.csv')}")
    print(f"Log: {os.path.join(output_dir, 'cluster.txt')}")
    print("=" * 60)


if __name__ == "__main__":
    main()
