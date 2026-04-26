#!/usr/bin/env bash
# AutoDock Vina redocking protocol for SN-38 vs Top1-DNA (PDB 1K4T)
# Run under an environment with: autodock-vina, openbabel, adfr-suite (for receptor prep)
set -eu
cd "$(dirname "$0")"

# 1) Split receptor / ligand
obabel 1K4T.pdb -O receptor.pdb --delete "HETATM and not resname TTG"
obabel 1K4T.pdb -O ligand.pdb   --select  "resname TTG"

# 2) Convert to .pdbqt
prepare_receptor -r receptor.pdb -o receptor.pdbqt
prepare_ligand   -l ligand.pdb   -o ligand.pdbqt

# 3) Define binding box around crystal ligand centroid
python3 -c "
from Bio.PDB import PDBParser
import numpy as np
s = PDBParser(QUIET=True).get_structure('l','ligand.pdb')
coords = np.array([a.coord for a in s.get_atoms()])
center = coords.mean(axis=0); dims = coords.ptp(axis=0) + 8
print(f'center_x={center[0]:.2f} center_y={center[1]:.2f} center_z={center[2]:.2f}')
print(f'size_x={dims[0]:.2f} size_y={dims[1]:.2f} size_z={dims[2]:.2f}')
" > box.conf

# 4) Run Vina
vina --receptor receptor.pdbqt --ligand ligand.pdbqt --config box.conf \
     --out sn38_vina_out.pdbqt --log vina.log --exhaustiveness 32 --num_modes 20

# 5) Extract top pose + compute RMSD vs crystal
obabel sn38_vina_out.pdbqt -O sn38_redocked.pdb
python3 -c "
from Bio.PDB import PDBParser; import numpy as np
cryst = np.array([a.coord for a in PDBParser(QUIET=True).get_structure('c','ligand.pdb').get_atoms()])
redoc = np.array([a.coord for a in PDBParser(QUIET=True).get_structure('r','sn38_redocked.pdb').get_atoms()])
n = min(len(cryst), len(redoc))
print(f'RMSD: {np.sqrt(((cryst[:n]-redoc[:n])**2).sum(axis=1).mean()):.3f} A')
"
