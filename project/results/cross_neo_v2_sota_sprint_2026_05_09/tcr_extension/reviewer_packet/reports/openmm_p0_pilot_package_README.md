# P0 OpenMM 10 ns Pilot Package

This package is for CUDA/OpenMM pilot MD of the two P0 TCR-pMHC complexes.

Recommended production pilot:

```bash
bash runpod_setup_openmm.sh
bash submit_p0_10ns_pilots.sh
```

Local CPU smoke test only:

```bash
python run_openmm_pilot.py --input inputs/<pilot>.minimized.pdb --outdir smoke/<pilot> --mode vacuum_smoke --platform CPU --steps 200 --report-steps 20 --minimize-iterations 200 --timestep-fs 0.5 --temperature-k 50
```

Claim boundary: 10 ns pilot MD is a structural stability diagnostic, not immunogenicity validation.
