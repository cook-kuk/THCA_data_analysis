"""ai_topvenue_mcp — top AI / stats / CV / NLP venue filter."""
from __future__ import annotations

TOP_AI_VENUES = frozenset({
    "NeurIPS", "ICML", "ICLR", "AAAI", "IJCAI",
    "CVPR", "ICCV", "ECCV",
    "ACL", "EMNLP", "NAACL",
    "KDD", "WWW", "SIGIR",
    "COLT", "AISTATS", "UAI",
    "MICCAI",  # medical imaging
})


def is_top_ai_venue(venue: str) -> bool:
    return venue.strip().upper() in {v.upper() for v in TOP_AI_VENUES}
