from .sources import (
    OpenAlex, CrossRef, PubMed, EuropePMC, SemanticScholar,
    Arxiv, BioRxiv, ICite, AFND, ClinicalTrials, Unpaywall,
)
from .tasks import (
    verify_references, competitive_landscape, hla_frequencies,
    misattribution_check, clinical_landscape, full_text_discovery,
    write_summary,
)

EMAIL = "kukshomr@gmail.com"
USER_AGENT = f"THCA-lit-enrich/1.0 (mailto:{EMAIL})"
