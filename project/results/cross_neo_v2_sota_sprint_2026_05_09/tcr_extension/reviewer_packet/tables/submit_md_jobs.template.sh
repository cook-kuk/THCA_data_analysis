#!/usr/bin/env bash
set -euo pipefail

# Template only. Fill INPUT_PDB after TCR-pMHC or pMHC structures pass QC.
# Recommended scratch: /data/thca/_tmp/md

# MDJOB_0001: HMTEVVRHC HLA-A*02:01 P0_MD_TCR_pMHC
# Required inputs: full_TCR_alpha_beta;MHC_sequence_or_template;TCR-pMHC_PDB
# Replicas: 3 x 100 ns
# INPUT_PDB=structures/MDJOB_0001.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0001
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0002: GADGVGKSAL HLA-C*08:02 P0_MD_TCR_pMHC
# Required inputs: full_TCR_alpha_beta;MHC_sequence_or_template;TCR-pMHC_PDB
# Replicas: 3 x 100 ns
# INPUT_PDB=structures/MDJOB_0002.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0002
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0003: MAWSLGVLVALPFPL HLA-B*40:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0003.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0003
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0004: MAWSLGVLVALPFPL HLA-A*02:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0004.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0004
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0005: FSLVFLVYSVFKNNV HLA-A*11:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0005.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0005
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0006: MAWSLGVLVALPFPL HLA-B*18:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0006.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0006
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0007: MAWSLGVLVALPFPL HLA-A*25:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0007.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0007
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0008: MSSLAATTFHWKKCR HLA-B*07:02 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0008.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0008
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0009: MSSLAATTFHWKKCR HLA-A*68:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0009.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0009
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0010: MSSLAATTFHWKKCR HLA-B*35:03 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0010.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0010
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0011: FSLVFLVYSVFKNNV HLA-B*38:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0011.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0011
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0012: FSLVFLVYSVFKNNV HLA-B*52:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0012.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0012
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0013: MIPSAAGIISLLDED HLA-B*15:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0013.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0013
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0014: MAWSLGVLVALPFPL HLA-C*03:04 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0014.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0014
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0015: PVLSGCSPRCLLQG HLA-B*15:01 P1_MD_pMHC_bulge
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0015.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0015
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0016: MSFSHLFYL HLA-A*02:01 P1_MD_pMHC_uncertain_TCR
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0016.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0016
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0017: VLMMPFSIV HLA-A*02:01 P1_MD_pMHC_uncertain_TCR
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0017.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0017
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0018: SSYYEPWGY HLA-A*01:01 P1_MD_pMHC_uncertain_TCR
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0018.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0018
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0019: LLDIVAPK HLA-A*11:01 P1_MD_pMHC_uncertain_TCR
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0019.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0019
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.

# MDJOB_0020: MIFFSIDRY HLA-A*26:01 P1_MD_pMHC_uncertain_TCR
# Required inputs: peptide_HLA_PDB;MHC_sequence_or_template
# Replicas: 3 x 50 ns
# INPUT_PDB=structures/MDJOB_0020.pdb
# OUTDIR=/data/thca/md_runs/MDJOB_0020
# mkdir -p "$OUTDIR"
# gmx pdb2gmx -f "$INPUT_PDB" -o "$OUTDIR/processed.gro" -water tip3p
# gmx editconf -f "$OUTDIR/processed.gro" -o "$OUTDIR/boxed.gro" -c -d 1.0 -bt dodecahedron
# gmx solvate -cp "$OUTDIR/boxed.gro" -cs spc216.gro -o "$OUTDIR/solv.gro" -p "$OUTDIR/topol.top"
# Add ions, minimize, equilibrate NVT/NPT, then run production replicas.
