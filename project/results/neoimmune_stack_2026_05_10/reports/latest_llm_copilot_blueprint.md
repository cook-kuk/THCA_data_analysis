# Latest LLM Copilot Blueprint

## Role
Use a current frontier LLM as a **rationale auditor and report copilot**, not as a scoring feature.

## Allowed LLM jobs
- Summarize candidate evidence into wet-lab-ready rationales.
- Detect contradictions: high score but low expression, HLA LOH, high WT similarity, weak presentation, leakage risk.
- Generate reviewer attack/defense tables.
- Draft patient-specific top-N handoff notes from structured TSV only.
- Convert model disagreement into experimental questions.

## Forbidden LLM jobs
- Directly assign immunogenicity labels.
- Override leakage flags.
- Convert public predictor scores into clean-track features.
- Claim clinical efficacy.

## Runtime contract
Set `NEOIMMUNE_LLM_MODEL` and provider keys externally. Store prompts and outputs with candidate IDs, input hashes, and model name for auditability.
