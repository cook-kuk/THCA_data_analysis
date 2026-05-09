# CROSS-Neo-TCR Literature and Tool Audit

Date: 2026-05-09

## Executive Decision

The TCR-aware layer should be built as an optional expert and diagnostic branch, not as a replacement for the main CROSS-Neo pMHC neoantigen ranker. The local registry contains substantial public TCR-pMHC evidence, but paired TCR alpha/beta coverage remains limited relative to the full neoantigen task and many public TCR records are pathogen-derived. Claims must therefore be limited to paired-TCR subsets, case interpretation, wetlab prioritization, and structure-grounded diagnostics until strict external paired-TCR benchmarks support broader claims.

## Local Availability Snapshot

| Component | Local status | Notes |
|---|---:|---|
| ColabFold / AlphaFold-Multimer wrapper | present | `colabfold_batch` found in `.venv/bin`; useful for AF2/Multimer-style batch structure generation. |
| AlphaFold Python package | present | `import alphafold` succeeds. |
| AlphaFold3 package | absent | `import alphafold3` fails; model parameter access has separate terms. |
| Boltz / Boltz-2 | absent | `import boltz` fails; no local install found. |
| Chai-1 | absent | `import chai_lab` fails; no local install found. |
| tFold / tFold-TCR | absent | `import tfold` fails; no local install found. |
| TCRdock / TCRmodel2 | absent | no local install found by filename scan; usable as future install or external workflow. |

## Tool Matrix

| Method / tool | Input requirements | Paired TCR alpha/beta required? | Peptide-HLA required? | MHC class | Structure output | Code available | Pretrained weights available | High-throughput suitability | Safe claim status | Integration priority |
|---|---|---:|---:|---|---:|---:|---:|---:|---|---:|
| TCRdock | TSV with organism, MHC class, MHC allele/sequence context, peptide, TRA/TRB V/J genes and CDR3s; can also parse ternary PDBs for docking geometry. | Yes for intended TCR-pMHC modeling. | Yes. | I and II parsing/modeling modes are documented. | Yes, via TCR-specialized AlphaFold setup plus docking geometry parsing. | Yes, GitHub. | Uses AlphaFold model parameters / fine-tuned workflow assets rather than a simple packaged checkpoint. | Medium; good for prioritized batches, not full-registry brute force without GPUs. | Diagnostic/structure-feature only until validated on heldout paired TCR-pMHC data. | P0 for structure diagnostics on exact paired matches; P1 for broader jobs. |
| TCRmodel2 | TCR alpha/beta variable sequences, peptide, MHC; web server and GitHub code. | Yes. | Yes. | Primarily TCR-pMHC class I/II modeling depending input/server support; verify per job. | Yes. | Yes, GitHub. | Uses adapted AlphaFold2/Multimer pipeline and model assets; not locally installed. | Medium; web/server or local setup, not as cheap as sequence-only models. | Diagnostic/structure-feature only. | P1. |
| AlphaFold-Multimer / ColabFold | FASTA sequences for all chains; for TCR-pMHC requires MHC alpha, beta2M for class I, peptide, TCR alpha/beta chains or reconstructed variable domains. | Yes for TCR-aware complex; no for pMHC-only. | Peptide + MHC sequence required for pMHC. | Protein complexes in general; TCR-pMHC class I/II if chains supplied. | Yes, PDB/mmCIF plus confidence metrics depending wrapper. | Yes, AlphaFold v2 code. | Yes for AF2/Multimer subject to model/data setup. | Medium; local ColabFold makes it practical for prioritized batches. | Diagnostic only for recognition unless linked to paired TCR labels. | P0/P1 because local ColabFold is available. |
| AlphaFold3 | JSON describing proteins and other molecules; for TCR-pMHC, all chains/peptide need to be specified. | Yes for TCR-pMHC. | Yes. | General biomolecular complexes; supports proteins and more molecule types. | Yes. | Inference code available; weights require Google access/terms. | Conditionally available via Google request/server terms. | Medium if installed/authorized; server has usage and molecule limits. | Diagnostic only; do not claim AF3 recognition correctness without benchmark. | P2 locally, P1 if weights/server access is approved. |
| Boltz / Boltz-2 | YAML/directory inputs describing biomolecules; can predict structures, and Boltz-2 also affinity-like outputs for supported ligand contexts. | Yes for TCR-pMHC if modeling full complex. | Yes for pMHC/TCR-pMHC. | General biomolecular complexes, not TCR-specific. | Yes. | Yes. | Yes, code and weights reported under MIT license. | High for installed GPU runs, but TCR-pMHC validation needed. | Diagnostic only for TCR recognition; promising but not locally active. | P2. |
| Chai-1 | FASTA containing proteins and other supported entities; optional MSAs/templates/restraints. | Yes for TCR-pMHC. | Yes. | General biomolecular complexes, not TCR-specific. | Yes, multiple samples. | Yes. | Yes, auto-downloaded weights in package workflow. | Medium/high on suitable GPUs. | Diagnostic only until TCR-pMHC benchmarked locally. | P2. |
| tFold-TCR | Full-chain TCR and pMHC sequences; supports TCR-only and TCR-pMHC predictors. | Yes for TCR-pMHC mode. | Yes for pMHC mode. | Trained on STCRDab/TCR-pMHC data; intended for TCR complexes. | Yes. | Yes, TencentAI4S/tfold. | Yes, downloadable model files referenced by project/Zenodo. | High; MSA-free workflow is designed for speed. | Diagnostic/structure-feature only until independently benchmarked on our splits. | P1, P0 if installed for exact paired subset. |
| STAG / STAG-LLM and structure-based GNN predictors | Modeled or experimental 3D TCR-pMHC structures plus sequence/PLM features depending variant. | Yes. | Yes. | TCR-pMHC; often class I-heavy depending dataset. | Consumes structure; does not primarily generate structure. | STAG datasets/code availability varies by release; STAG-LLM paper is open. | Not established as a packaged local checkpoint here. | Low/medium; structure generation dominates cost. | Diagnostic/modeling research only; strong supplement candidate, not main claim today. | P2. |
| UniPMT | Peptide, MHC pseudo-sequence/HLA, TCR embeddings in a heterogeneous graph for P-M, P-T, and P-M-T tasks. | TCR required for P-T/P-M-T; paired alpha/beta support must be verified from model data format. | Yes for P-M-T, optional for P-T. | HLA/pMHC task dependent. | No. | Yes. | Model files referenced in repository. | High after setup; GPU recommended. | Baseline/expert only; claim only on strict paired-TCR subset. | P1/P2. |
| UnifyImmun | Peptide plus HLA and/or TCR sequences; cross-attention model for pHLA and pTCR binding. | TCR required for pTCR; paired alpha/beta detail depends task/config. | HLA required for pHLA; pTCR can be peptide+TCR. | pHLA and pTCR tasks; class support should follow training data. | No. | Yes. | Trained model folder in repository. | High after GPU setup. | Baseline/expert only; not evidence to infer TCR recognition from peptide-HLA alone. | P1/P2. |
| TCRBagger | Patient peptide list, TCR list or constructed bags; optional RNA-seq, VCF, HLA alleles; integrates netMHCpan/MuPeXI-style neoantigen workflow. | No; bagged TCR repertoire rather than paired cognate alpha/beta TCR-pMHC. | HLA-I peptide candidate context required for neoantigen application. | HLA-I. | No. | Yes. | Repository includes `Models` directory. | Medium; older TensorFlow/MiXCR/netMHCpan dependencies. | Wetlab prioritization/diagnostic only; weak labels are not cognate TCR ground truth. | P2. |
| TPepRet | TCR and peptide sequences for TCR-peptide binding characterization. | Likely sequence-level TCR input; paired alpha/beta support must be checked in code before use. | Peptide required; explicit HLA dependence is not central in title/abstract. | Not clearly HLA-stratified from audit source. | No. | Yes. | Code available; checkpoint packaging must be checked. | High after setup. | Benchmark baseline only; diagnostic if no HLA is modeled. | P2. |
| ERGO-II | TCR beta CDR3 and peptide always; optional TCR alpha, V/J genes, MHC, and CD4/CD8 type. | No; alpha optional and model has missing-alpha path. | Peptide required; MHC optional. | MHC optional, class via data/features. | No. | Yes. | Two models in repo for McPAS and VDJdb; web tool also available. | High; web tool capped at 50,000 pairs, local model usable for larger batches. | Good sequence baseline; claims only under strict TCR/epitope/study split. | P1. |
| NetTCR-2.0 | CSV with CDR3a, CDR3b, peptide, binder for alpha+beta model; beta-only training sets also provided. | Preferred yes for NetTCR-2.0 alpha+beta; beta-only models exist. | Peptide required; HLA not explicit in the basic file format. | Often peptide/TCR rather than explicit HLA. | No. | Yes. | Pretrained web server; repo has training data/code. | High. | P0/P1 baseline on paired-TCR subset; no peptide-HLA-only inference. | P1. |
| ImRex | TCR beta CDR3 and epitope sequence encoded as interaction map for CNN. | No; beta-chain focused. | Peptide/epitope required; HLA not central. | Not explicit. | No. | Yes. | Small pretrained models included. | High. | Historical sequence baseline; diagnostic only for neoantigen recognition. | P2. |
| pMTnet | CSV with CDR3 beta, antigen peptide, HLA; trained model/library with HLA sequences and background TCRs. | No; beta CDR3 only. | Yes: peptide + HLA. | HLA class context from library, largely class I-oriented in neoantigen use. | No. | Yes. | Library/trained models provided in release/workflow. | High after legacy dependency setup. | Useful diagnostic feature; do not treat beta-only prediction as paired cognate recognition. | P1/P2. |
| TITAN | TCR sequence table, epitope table, train/test label pairs; epitope can be amino acid or SMILES-style representation. | No; TCR sequence input, generally beta-chain style datasets. | Peptide/epitope required; HLA not explicit. | Not explicit. | No. | Yes. | `trained_model` directory available. | High after setup. | Historical benchmark baseline; diagnostic only. | P2. |

## Database / Dataset Matrix

| Database / dataset | Input / contents | Paired TCR alpha/beta required? | Peptide-HLA required? | MHC class | Structure output | Code / access | Pretrained weights | High-throughput suitability | Safe claim status | Integration priority |
|---|---|---:|---:|---|---:|---:|---:|---:|---|---:|
| VDJdb | Curated antigen-specific TCR records with CDR3, V/J genes, peptide, antigen source, MHC, method metadata. | No; contains alpha-only, beta-only, and paired records. | Often yes, but missing/partial records exist. | I and II plus non-human MHC contexts. | Some records link structures. | Public database/export. | N/A. | High. | P0 evidence source; diagnostic unless label/split compatibility is proven. | P0. |
| McPAS-TCR | Manually curated pathology-associated TCR sequences with disease/pathology, antigen/epitope when available. | No; many beta-only records. | No; many pathology-associated rows lack full peptide-HLA. | Mixed/partial. | No primary structure output. | Public database/export. | N/A. | High. | Diagnostic/public-overlap evidence; pathology association is not cognate neoantigen recognition. | P0. |
| IEDB | T-cell, B-cell, and MHC ligand/epitope assays; receptor sequence exports and IEDB-3D structures exist for some entries. | No. | Often peptide-HLA for T-cell/MHC assays; not universal. | I and II. | IEDB-3D / PDB links for subset. | Public database/API/export. | N/A. | High. | P0 peptide-HLA/T-cell assay evidence; labels not transferred without assay/source compatibility. | P0. |
| PIRD | Pan immune repertoire repository with raw/processed TCR/BCR sequences across species/phenotypes. | No; paired chains depend on study technology. | No. | Not inherently pMHC-restricted. | No primary structure output. | Public database. | N/A. | High for repertoire mining; lower for cognate TCR-pMHC. | Registry expansion only; not recognition ground truth. | P3 until local source is acquired. |
| 10x Genomics immune profiling datasets | Single-cell V(D)J, gene expression, surface/antigen feature barcode datasets; selected public datasets link paired alpha/beta TCRs to pMHC multimers/dextramers. | Often yes when V(D)J + antigen capture used. | Antigen-capture panels provide peptide-HLA multimer labels; not all datasets. | Usually class I panels in public antigen examples, but technology is broader. | No direct structure output. | Public support datasets / publications. | N/A. | High once downloaded and QC-filtered. | P0 for paired-TCR antigen labels, but multimer background must be filtered. | P0/P1. |
| TetTCR-seq / TetTCR-SeqHD datasets | DNA-barcoded pMHC tetramer binding linked to single-cell TCR sequence and sometimes expression/surface phenotype. | Yes by design for single-cell TCR recovery, subject to QC. | Yes, tetramer pMHC identity. | Mostly class I CD8 workflows in cited studies. | No direct structure output. | Publication datasets. | N/A. | Medium/high; excellent for cognate specificity after QC. | Strong TCR-aware subset evidence; still benchmark/diagnostic unless external heldout supports claims. | P0/P1. |
| STCRDab | Structural TCR database with TCR structures and annotations including docking angles. | Yes for alpha/beta TCR structures where available. | Some entries are TCR-pMHC; not all TCR-only structures require peptide-HLA. | I, II, and other immune recognition contexts. | Yes, experimental PDB structures/annotations. | Public database. | N/A. | Medium; finite structure set. | Structure benchmarking and feature calibration; not broad prediction ground truth. | P0 for structure QC, P1 for modeling templates. |
| TCR3d | Structural repertoire database of TCR-antigen complexes, including class I, class II, CD1d, MR1, and gamma-delta contexts. | Yes for TCR complex structures. | For pMHC complexes yes; also includes MHC-like antigen contexts. | I, II, CD1d, MR1, gamma-delta. | Yes, experimentally determined structures and geometry annotations. | Public database. | N/A. | Medium; curated structural set. | P0/P1 structure evidence and geometry validation, not high-throughput label source. | P0/P1. |

## Practical CROSS-Neo-TCR Use

| Use case | Recommended resources | Claim boundary |
|---|---|---|
| TCR-available subset training | VDJdb, 10x antigen capture, TetTCR-seq, IEDB receptor exports, selected McPAS rows | Only paired alpha/beta rows with compatible peptide-HLA labels should enter Tier 1 modeling. |
| Baseline sequence models | NetTCR-2.0, ERGO-II, pMTnet, TITAN, ImRex, TPepRet, UniPMT, UnifyImmun | Baselines for strict group splits; not universal neoantigen recognition predictors. |
| Structure feature generation | TCRdock, AlphaFold-Multimer/ColabFold, tFold-TCR, TCRmodel2, AF3 if authorized, Boltz/Chai if installed | Use confidence/contact/geometry features for diagnostics and wetlab prioritization. |
| Structure benchmarking/calibration | STCRDab, TCR3d, TCRdock parser, STAG/STAG-LLM/TCRLens literature | Benchmark and explain cases; avoid main SOTA claims without clean external tests. |
| Repertoire-level patient prioritization | TCRBagger, 10x repertoire data, local tumor TCR-seq if available | Patient/repertoire diagnostic only; does not identify cognate TCR-pMHC pairs by itself. |

## Immediate Integration Decisions

- P0 data sources for the registry: VDJdb, McPAS-TCR, IEDB, local 10x/VDJdb chunks, NEPdb-derived neoantigen/TCR fields.
- P0/P1 linkage evidence: exact peptide+HLA+paired TCR is strongest; peptide-only and near-peptide matches remain non-definitive.
- P1 modeling baselines: ERGO-II, NetTCR-2.0, pMTnet, UniPMT/UnifyImmun after environment setup.
- P1 structure workflow: use local ColabFold/AlphaFold first; add tFold-TCR or TCRdock if installation time is justified by paired exact matches and wetlab candidates.
- P2/P3 future work: Boltz-2, Chai-1, AlphaFold3 local inference, STAG-style GNNs, PIRD expansion.

## Sources

- TCRdock GitHub and documentation: https://github.com/phbradley/TCRdock
- TCRdock paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9859041/
- TCRmodel2 GitHub: https://github.com/piercelab/tcrmodel2
- TCRmodel2 paper: https://academic.oup.com/nar/article/51/W1/W569/7151345
- AlphaFold v2 / Multimer GitHub: https://github.com/google-deepmind/alphafold
- AlphaFold3 GitHub: https://github.com/google-deepmind/alphafold3
- AlphaFold3 paper: https://www.nature.com/articles/s41586-024-07487-w
- Boltz GitHub: https://github.com/jwohlwend/boltz
- Chai-1 GitHub: https://github.com/chaidiscovery/chai-lab
- tFold / tFold-TCR GitHub: https://github.com/TencentAI4S/tfold
- STAG-LLM paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC12800378/
- STAG publication page: https://kavrakilab.rice.edu/publications/slone2025-tcr-binding.html
- TCRLens structure-aware benchmark: https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbag066/8496266
- UniPMT GitHub: https://github.com/ethanmock/UniPMT
- UnifyImmun GitHub: https://github.com/hliulab/UnifyImmun
- UnifyImmun paper: https://www.nature.com/articles/s42256-024-00973-w
- TCRBagger GitHub: https://github.com/bm2-lab/TCRBagger
- TPepRet paper/code link: https://academic.oup.com/bioinformatics/article/41/1/btaf022/7989444
- ERGO-II GitHub: https://github.com/IdoSpringer/ERGO-II
- NetTCR-2.0 GitHub: https://github.com/mnielLab/NetTCR-2.0
- ImRex GitHub: https://github.com/pmoris/ImRex
- pMTnet GitHub: https://github.com/tianshilu/pMTnet
- TITAN GitHub: https://github.com/PaccMann/TITAN
- TCR-pMHC model comparison: https://pmc.ncbi.nlm.nih.gov/articles/PMC10152969/
- VDJdb paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC5753233/
- VDJdb 2019 update: https://pmc.ncbi.nlm.nih.gov/articles/PMC6943061/
- McPAS-TCR paper: https://academic.oup.com/bioinformatics/article/33/18/2924/3803440
- IEDB receptor sequence export docs: https://help.iedb.org/hc/en-us/articles/115002891972-TCR-and-Antibody-Sequence-Data
- IEDB epitope scope docs: https://help.iedb.org/hc/en-us/articles/114094147471-IEDB-Epitopes
- PIRD paper: https://academic.oup.com/bioinformatics/article/36/3/897/5543102
- PIRD database page: https://db.cngb.org/pird/
- 10x paired immune repertoire overview: https://www.10xgenomics.com/blog/high-definition-immunology-paired-immune-repertoire-profiling-at-single-cell-resolution
- 10x antigen-recognition overview: https://www.10xgenomics.com/blog/a-sequencing-approach-to-t-cell-receptor-antigen-recognition
- TetTCR-seq paper: https://pubmed.ncbi.nlm.nih.gov/30418433/
- TetTCR-SeqHD paper: https://www.nature.com/articles/s41590-021-01073-2
- STCRDab paper: https://academic.oup.com/nar/article/46/D1/D406/4566020
- TCR3d paper: https://academic.oup.com/bioinformatics/article/35/24/5323/5523179
