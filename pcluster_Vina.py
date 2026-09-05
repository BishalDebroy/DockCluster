"""
dockcluster3d.py

A production-quality script for 3D spatial clustering of ligand PDBQT files
based on center-of-mass (COM) coordinates.

It reads all .pdbqt files from an input directory, computes the mass-weighted
COM of each ligand, clusters them using a connectivity-based algorithm
(transitive closure with a distance cutoff), and moves each ligand into a
cluster subdirectory. An optional target-coordinate mode reorders the
clusters so that the cluster closest to a specific point becomes cluster_1.

Usage examples:
    # Normal clustering
    python dockcluster3d.py --input ./pdbqt_files --cutoff 1.0

    # Target-coordinate mode
    python dockcluster3d.py --input ./pdbqt_files --target 12.5 8.9 -3.1 --cutoff 1.0

    # Show the script version
    python dockcluster3d.py --version

Dependencies: None (pure Python 3.6+ standard library)

Version: 1.0

Changelog:
    1.0 - Clusters are now ordered by decreasing number of ligands (largest
          cluster first, smallest last) in both cluster_log.txt and the
          cluster_N output directories. In target-coordinate mode, cluster_1
          is still always the cluster closest to the target; the remaining
          clusters are ordered by decreasing size after that.
"""

import os
import sys
import argparse
import math
import shutil
import glob
import warnings
from collections import defaultdict
from typing import List, Tuple, Dict, Optional, NamedTuple

__version__ = "1.0"

# ---------------------------
# Atomic masses (in g/mol)
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
    "Na": 22.990,
    "Mg": 24.305,
    "K": 39.098,
    "Ca": 40.078,
    "Fe": 55.845,
    "Zn": 65.38,
    "Se": 78.971,
    "Mn": 54.938,
}

# Case-insensitive lookup: "cl" / "CL" / "Cl" all resolve to the "Cl" key.
_ATOMIC_MASS_CI: Dict[str, str] = {sym.lower(): sym for sym in ATOMIC_MASS}


def _lookup_element(candidate: str) -> Optional[str]:
    """
    Case-insensitively resolve a candidate string to a known element symbol
    from ATOMIC_MASS. Returns the canonically-cased symbol, or None if the
    candidate isn't a recognized element.
    """
    if not candidate:
        return None
    return _ATOMIC_MASS_CI.get(candidate.lower())


class Atom(NamedTuple):
    """Representation of a single atom from a PDBQT file."""
    x: float
    y: float
    z: float
    element: str
    mass: float


class Ligand(NamedTuple):
    """Information about one ligand."""
    filename: str
    com: Tuple[float, float, float]  # (x, y, z)
    atoms: List[Atom]


# -----------------------------
# PDBQT parsing
def parse_pdbqt(filepath: str) -> List[Atom]:
    """
    Parse a PDBQT file and return a list of Atom objects.

    Reads ATOM and HETATM records, extracts coordinates using fixed-column
    positions (following the standard PDB/AutoDock-PDBQT layout), and
    determines the element symbol. If the element cannot be determined, a
    warning is issued and the atom is skipped.
    """
    atoms = []
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue

            try:
                # Standard PDB columns (1-indexed): x=31-38, y=39-46, z=47-54
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                warnings.warn(f"Invalid coordinates in {filepath}, skipping line.")
                continue

            # Determine element symbol
            element = None

            # 1) Try the AutoDock atom-type / element column (columns 78-79)
            elem_col = line[77:79].strip()
            if elem_col:
                candidate = elem_col.strip('0123456789+-')
                element = _lookup_element(candidate)

            # 2) If not found, infer from the atom name (columns 13-16)
            if element is None:
                atom_name = line[12:16].strip()
                letters = ''.join(ch for ch in atom_name if ch.isalpha())
                if letters:
                    # Prefer the single-letter interpretation first, since
                    # protein atom names like "CA" (alpha-carbon) should
                    # resolve to Carbon, not Calcium.
                    element = _lookup_element(letters[0])
                    if element is None and len(letters) >= 2:
                        element = _lookup_element(letters[:2])

            # 3) If still None, the residue name sometimes contains the element
            if element is None:
                res_name = line[17:20].strip()
                element = _lookup_element(res_name)

            # 4) If still no element, skip this atom with a warning.
            if element is None:
                warnings.warn(f"Could not determine element for atom in {filepath}. Skipping atom.")
                continue

            mass = ATOMIC_MASS.get(element)
            if mass is None:
                warnings.warn(f"Unknown element '{element}' in {filepath}. Skipping atom.")
                continue

            atoms.append(Atom(x=x, y=y, z=z, element=element, mass=mass))

    if not atoms:
        warnings.warn(f"No atoms found in {filepath}.")
    return atoms


# -----------------------------
# Center of mass calculation
def compute_com(atoms: List[Atom]) -> Tuple[float, float, float]:
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


def distance_3d(p1: Tuple[float, float, float],
                 p2: Tuple[float, float, float]) -> float:
    return math.sqrt((p1[0] - p2[0]) ** 2 +
                      (p1[1] - p2[1]) ** 2 +
                      (p1[2] - p2[2]) ** 2)


# -----------------------------
# Connectivity-based clustering
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


def cluster_ligands(ligands: List[Ligand], cutoff: float) -> List[List[int]]:
    """
    Cluster ligands based on COM distances using transitive closure.
    Returns a list of clusters (each a list of ligand indices).
    """
    n = len(ligands)
    uf = UnionFind(n)

    # Compute pairwise distances and union if <= cutoff
    for i in range(n):
        com_i = ligands[i].com
        for j in range(i + 1, n):
            com_j = ligands[j].com
            if distance_3d(com_i, com_j) <= cutoff:
                uf.union(i, j)

    # Group indices by root
    clusters_dict: Dict[int, List[int]] = defaultdict(list)
    for idx in range(n):
        root = uf.find(idx)
        clusters_dict[root].append(idx)

    return list(clusters_dict.values())


# -----------------------------
# Determining ordering of clusters
def cluster_center_of_mass(ligands: List[Ligand], indices: List[int]) -> Tuple[float, float, float]:
    """Compute the mean of the ligand COMs in the cluster."""
    if not indices:
        return (0.0, 0.0, 0.0)
    sum_x = sum_y = sum_z = 0.0
    for idx in indices:
        x, y, z = ligands[idx].com
        sum_x += x
        sum_y += y
        sum_z += z
    n = len(indices)
    return (sum_x / n, sum_y / n, sum_z / n)


def order_clusters(clusters: List[List[int]],
                    ligands: List[Ligand]) -> List[List[int]]:
    """
    Sort clusters by decreasing number of ligands (largest cluster first,
    smallest last). Ties (equal-sized clusters) are broken deterministically
    by cluster COM (X, then Y, then Z).
    """
    cluster_coms = []
    for cluster in clusters:
        com = cluster_center_of_mass(ligands, cluster)
        cluster_coms.append(com)
    # Sort by cluster size descending; break ties by (com_x, com_y, com_z)
    sorted_pairs = sorted(zip(cluster_coms, clusters),
                           key=lambda pair: (-len(pair[1]), pair[0][0], pair[0][1], pair[0][2]))
    return [pair[1] for pair in sorted_pairs]


# -----------------------------
# Target mode: find cluster closest to target
def find_closest_cluster_to_target(ligands: List[Ligand],
                                    clusters: List[List[int]],
                                    target: Tuple[float, float, float],
                                    radius: float) -> Tuple[Optional[int], Optional[List[Dict]]]:
    """
    Identify the cluster that contains a ligand within radius of the target.
    """
    # Gather all matching ligands: ligand index, distance, cluster index
    matches = []  # (ligand_idx, distance, cluster_idx)
    for c_idx, cluster in enumerate(clusters):
        for lig_idx in cluster:
            dist = distance_3d(ligands[lig_idx].com, target)
            if dist <= radius:
                matches.append((lig_idx, dist, c_idx))

    if not matches:
        return None, None

    # Choose the match with smallest distance
    matches.sort(key=lambda x: x[1])  # sort by distance
    best_lig_idx, best_dist, best_c_idx = matches[0]

    # Build info for logging (all matching molecules)
    match_info = []
    for lig_idx, dist, c_idx in matches:
        match_info.append({
            'ligand_idx': lig_idx,
            'filename': ligands[lig_idx].filename,
            'distance': dist,
            'cluster_idx': c_idx,
            'is_closest': (lig_idx == best_lig_idx and c_idx == best_c_idx),
        })

    return best_c_idx, match_info


# -----------------------------
# Moving files and writing log
def move_ligands_to_clusters(ligands: List[Ligand],
                              ordered_clusters: List[List[int]],
                              output_dir: str,
                              cluster_prefix: str = "cluster_") -> None:
    """Move ligand files into numbered cluster subdirectories."""
    for i, cluster in enumerate(ordered_clusters, start=1):
        cluster_dir = os.path.join(output_dir, f"{cluster_prefix}{i}")
        os.makedirs(cluster_dir, exist_ok=True)
        for idx in cluster:
            src = ligands[idx].filename
            # src may be a full path; use the base name for the destination
            base = os.path.basename(src)
            dst = os.path.join(cluster_dir, base)
            shutil.move(src, dst)


def write_cluster_log(ligands: List[Ligand],
                       ordered_clusters: List[List[int]],
                       output_dir: str,
                       cutoff: float,
                       input_dir: str,
                       target: Optional[Tuple[float, float, float]] = None,
                       radius: Optional[float] = None,
                       target_cluster_idx: Optional[int] = None,
                       match_info: Optional[List[Dict]] = None,
                       skipped_files: Optional[List[str]] = None) -> None:
    """Write cluster_log.txt with full details."""
    log_path = os.path.join(output_dir, "cluster_log.txt")
    with open(log_path, 'w') as f:

        # Header summary
        f.write("=" * 60 + "\n")
        f.write("LIGAND COM CLUSTERING REPORT\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Script version         : {__version__}\n")
        f.write(f"Input directory        : {input_dir}\n")
        f.write(f"Output directory       : {output_dir}\n")
        f.write(f"Clustering cutoff      : {cutoff:.3f} \u00c5\n")
        f.write(f"Number of PDBQT files  : {len(ligands)}\n")
        f.write(f"Files skipped          : {len(skipped_files) if skipped_files else 0}\n")
        f.write(f"Number of clusters     : {len(ordered_clusters)}\n")
        if ordered_clusters:
            sizes = [len(c) for c in ordered_clusters]
            f.write(f"Largest cluster        : {max(sizes)} ligands\n")
            f.write(f"Smallest cluster       : {min(sizes)} ligands\n")
        else:
            f.write("Largest cluster        : 0\n")
            f.write("Smallest cluster       : 0\n")

        # If target mode, add target info
        if target is not None:
            f.write("\n" + "=" * 60 + "\n")
            f.write("TARGET COORDINATE MODE\n")
            f.write("=" * 60 + "\n")
            f.write(f"Target coordinates     : {target[0]:.3f}, {target[1]:.3f}, {target[2]:.3f} \u00c5\n")
            f.write(f"Target radius          : {radius:.3f} \u00c5\n")
            if target_cluster_idx is not None:
                # The chosen number after renumbering
                f.write("\nSelected cluster       : cluster_1\n")
                if match_info:
                    closest = [m for m in match_info if m.get('is_closest')]
                    if closest:
                        m = closest[0]
                        f.write(f"Reason                 : {m['filename']} COM is {m['distance']:.3f} \u00c5 from target.\n")
                        # List all matching ligands
                        f.write("\nAll matching ligands within radius:\n")
                        # Built as a variable rather than inline: a nested
                        # string literal containing a "\u00c5" escape inside
                        # an f-string's {...} expression is a SyntaxError on
                        # Python < 3.12 (PEP 701 lifted this restriction).
                        distance_header = "Distance (\u00c5)"
                        f.write(f"{'Filename':<35}{distance_header:>15}\n")
                        f.write("=" * 50 + "\n")
                        for m in sorted(match_info, key=lambda x: x['distance']):
                            f.write(f"{m['filename']:<35}{m['distance']:>15.3f}\n")
                        # Indicate that multiple clusters may have had matches
                        f.write("\nNote: The cluster chosen as cluster_1 contains the ligand with the smallest distance.\n")
                        unique_clusters = set(m['cluster_idx'] for m in match_info)
                        if len(unique_clusters) > 1:
                            f.write("WARNING: Multiple clusters contain ligands within the radius.\n")
                            f.write("Only the cluster with the closest ligand was selected as cluster_1.\n")
                            f.write("All other clusters were renumbered starting from cluster_2.\n")
            else:
                f.write("\nNo ligand found within the target radius.\n")
                f.write("Clusters numbered normally (no reordering).\n")
            f.write("\n")

        # Now list each cluster
        for i, cluster in enumerate(ordered_clusters, start=1):
            com = cluster_center_of_mass(ligands, cluster)
            f.write("\n" + "=" * 60 + "\n")
            f.write(f"CLUSTER {i}\n")
            f.write(f"Number of ligands: {len(cluster)}\n")
            f.write(f"Cluster COM   : X = {com[0]:.3f} \u00c5, Y = {com[1]:.3f} \u00c5, Z = {com[2]:.3f} \u00c5\n")
            f.write('\n')
            f.write("-" * 80 + "\n")
            f.write(f"{'Filename':<35}{'X_COM':>10}{'Y_COM':>10}{'Z_COM':>10}{'Distance to cluster COM':>25}\n")
            f.write("-" * 80 + "\n")
            for idx in cluster:
                lig = ligands[idx]
                com_lig = lig.com
                dist_to_cluster = distance_3d(com_lig, com)
                f.write(f"{os.path.basename(lig.filename):<35}"
                        f"{com_lig[0]:>10.3f}"
                        f"{com_lig[1]:>10.3f}"
                        f"{com_lig[2]:>10.3f}"
                        f"{dist_to_cluster:>25.3f}\n")
            f.write("\n")

        # If any files were skipped, list them
        if skipped_files:
            f.write("\n" + "=" * 60 + "\n")
            f.write("SKIPPED FILES\n")
            f.write("=" * 60 + "\n")
            for sf in skipped_files:
                f.write(f"  {sf}\n")
            f.write("\n")


# -----------------------------
# Main CLI
def main():
    parser = argparse.ArgumentParser(
        description="Cluster PDBQT ligands based on center-of-mass distance.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--version", action="version",
                         version=f"%(prog)s {__version__}",
                         help="Show the script version and exit")
    parser.add_argument("--input", required=True, help="Directory containing .pdbqt files")
    parser.add_argument("--output", default=None, help="Output directory (default: <input>/clustered)")
    parser.add_argument("--cutoff", type=float, default=1.0, help="Clustering cutoff in \u00c5 (default: 1.0)")
    parser.add_argument("--target", nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                         help="Target coordinate (X Y Z) for target-coordinate mode")
    parser.add_argument("--radius", type=float, default=1.0, help="Radius around target for selection (default: 1.0)")

    args = parser.parse_args()

    input_dir = args.input
    if not os.path.isdir(input_dir):
        print(f"ERROR: input directory '{input_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    if args.output is None:
        output_dir = os.path.join(input_dir, "clustered")
    else:
        output_dir = args.output

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    cutoff = args.cutoff
    if cutoff < 0:
        print("ERROR: Cutoff must be non-negative.", file=sys.stderr)
        sys.exit(1)

    target = None
    radius = None
    if args.target is not None:
        target = tuple(args.target)
        radius = args.radius
        if radius < 0:
            print("ERROR: Radius must be non-negative.", file=sys.stderr)
            sys.exit(1)

    # Gather PDBQT files. glob is case-sensitive on Linux/macOS, so also
    # pick up mixed-case extensions like .PDBQT, deduplicating the result.
    files = sorted(set(
        glob.glob(os.path.join(input_dir, "*.pdbqt")) +
        glob.glob(os.path.join(input_dir, "*.PDBQT")) +
        glob.glob(os.path.join(input_dir, "*.Pdbqt"))
    ))
    if not files:
        print(f"Warning: No .pdbqt files found in '{input_dir}'.", file=sys.stderr)
        sys.exit(0)

    # Process each file
    ligands = []
    skipped = []
    for fpath in files:
        try:
            atoms = parse_pdbqt(fpath)
            if not atoms:
                skipped.append(fpath)
                continue
            com = compute_com(atoms)
            ligands.append(Ligand(filename=fpath, com=com, atoms=atoms))
        except Exception as e:
            warnings.warn(f"Error processing {fpath}: {e}")
            skipped.append(fpath)

    if not ligands:
        print("No valid ligands found. Exiting.", file=sys.stderr)
        with open(os.path.join(output_dir, "cluster_log.txt"), 'w') as f:
            f.write("No valid ligands processed.\n")
        sys.exit(1)

    # Perform clustering
    raw_clusters = cluster_ligands(ligands, cutoff)
    ordered_clusters = order_clusters(raw_clusters, ligands)

    # Target mode handling
    target_cluster_idx = None
    match_info = None
    if target is not None:
        best_c_idx, match_info = find_closest_cluster_to_target(ligands, ordered_clusters, target, radius)
        if best_c_idx is not None:
            selected_cluster = ordered_clusters.pop(best_c_idx)
            ordered_clusters.insert(0, selected_cluster)
            target_cluster_idx = 0
        else:
            target_cluster_idx = None

    # Move ligand files into cluster directories
    move_ligands_to_clusters(ligands, ordered_clusters, output_dir)

    # Write log
    write_cluster_log(
        ligands=ligands,
        ordered_clusters=ordered_clusters,
        output_dir=output_dir,
        cutoff=cutoff,
        input_dir=input_dir,
        target=target,
        radius=radius if target is not None else None,
        target_cluster_idx=target_cluster_idx,
        match_info=match_info,
        skipped_files=skipped
    )

    # Console summary
    print("=" * 60)
    print("Ligand COM Clustering Complete")
    print("=" * 60)
    print(f"Script version           : {__version__}")
    print(f"Input files processed    : {len(ligands)}")
    print(f"Files skipped            : {len(skipped)}")
    print(f"Clusters generated       : {len(ordered_clusters)}")
    print(f"Clustering cutoff        : {cutoff:.3f} \u00c5")
    print(f"Output directory         : {output_dir}")
    print(f"Log file                 : {os.path.join(output_dir, 'cluster_log.txt')}")
    print("=" * 60)


if __name__ == "__main__":
    main()
