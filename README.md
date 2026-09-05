# DockCluster: COM‑based spatial clustering of docking poses.

DockCluster provides two standalone Python scripts that perform 3D clustering of ligand conformations based **only** on the **mass‑weighted center of mass (COM)** of each ligand. Clustering is transitive (connectivity‑based), with a default cutoff of **1.0 Å**.

## Features

- **For Vina (PDBQT)**:  
  - Reads all `.pdbqt` files from a directory.  
  - Calculates the mass‑weighted COM of each ligand.  
  - Clusters by COM‑to‑COM Euclidean distance (transitive).  
  - Moves each ligand into a numbered `cluster_N/` subdirectory.  
  - Generates a detailed `cluster_log.txt` with cluster COM and per‑ligand coordinates.

- **For AutoDock4 (DLG)**:  
  - Parses a `.dlg` docking log file.  
  - Extracts each docking conformation’s coordinates and the six energy terms.  
  - Clusters the conformations by COM‑to‑COM distance.  
  - For each cluster, selects the representative conformation with the **lowest Estimated Free Energy of Binding**.  
  - Produces `cluster.csv` (one row per cluster, with all six energies from the representative) and `cluster.txt` (full log).

## Installation

No external libraries are required. Ensure you have **Python 3.6 or later**.

```bash
# Clone the repository
git clone https://github.com/yourusername/DockCluster3D.git
cd DockCluster3D
```

## Clustering of Vina outputs
1. If target binding site is unknown,
```
python pcluster_Vina.py --input ./pdbqt_files --cutoff 2.0
```
2. if target binding site is known (eg. x=3.21 y=-5.24 z=17.92),
```
python pcluster_Vina.py --input ./pdbqt_files --target 3.21 5.24 17.92 --radius 1.0
```
3. For more options,
```
python pcluster_Vina.-py --help
```

## Clustering of AD4 outputs
1. If target binding site is unknown
```
python pcluster_AD4.py --input docking.dlg --cutoff 2.0
```

If you use DockCluster in your research, please cite:
```
@misc{BishalDockCluster,
  author       = {Bishal},
  title        = {DockCluster3D: COM‑based spatial clustering of docking poses},
  year         = {2026},
  howpublished = {\url{https://github.com/BishalDebroy/Dockcluster}},
  note         = {MIT License}
}
```
