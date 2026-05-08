"""CLI entry point for Cancer Discovery OS.

Usage examples:

    python3 -m agentic_research.domains.cancer_discovery_os.cli manifest
    python3 -m agentic_research.domains.cancer_discovery_os.cli search PTC
    python3 -m agentic_research.domains.cancer_discovery_os.cli ask "plan a CRISPR screen for thyroid dedifferentiation"
    python3 -m agentic_research.domains.cancer_discovery_os.cli export --format markdown > lab_assets.md
    python3 -m agentic_research.domains.cancer_discovery_os.cli export --format json > lab_assets.json

The export subcommand is what you attach to a grant proposal as evidence
of platform + asset readiness.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from enum import Enum

from .manifest import EntryKind, default_yu_cook_manifest
from .runner import CancerDiscoveryOS, OSQuery


def _entry_to_dict(e):
    d = asdict(e)
    d["kind"] = e.kind.value
    return d


def _cmd_manifest(args) -> int:
    m = default_yu_cook_manifest()
    print(f"Lab manifest — PI={m.pi}, Co-PI={m.co_pi}, entries={len(m)}")
    by_kind: dict[str, int] = {}
    for e in m.entries:
        by_kind[e.kind.value] = by_kind.get(e.kind.value, 0) + 1
    for k, n in sorted(by_kind.items()):
        print(f"  {k:11s} {n}")
    if args.verbose:
        print()
        for e in m.entries:
            tag = f"[{e.kind.value}/{e.status}]"
            print(f"  {tag:22s} {e.id:32s} {e.title}")
    return 0


def _cmd_search(args) -> int:
    os_ = CancerDiscoveryOS()
    r = os_.handle(OSQuery(f"asset lookup {args.term}", payload={"search": args.term}))
    print(f"intent={r.intent.value} hits={len(r.result)}")
    for e in r.result:
        print(f"  [{e.kind.value}] {e.id} — {e.title}")
        if args.verbose and e.summary:
            print(f"      {e.summary}")
    return 0


def _cmd_ask(args) -> int:
    os_ = CancerDiscoveryOS()
    r = os_.handle(OSQuery(args.text))
    print(f"intent={r.intent.value} handler={r.handler}")
    if r.notes:
        for n in r.notes:
            print(f"  note: {n}")
    if r.result is not None:
        if hasattr(r.result, "ranking"):
            print(f"  topics={len(r.result.topics)} battles={len(r.result.battle_log)}")
            for t in r.result.ranking[:5]:
                sc = r.result.final_scores[t.title]
                print(f"  [{sc.tier} {sc.total}] {t.title}")
        elif isinstance(r.result, list):
            for e in r.result[:10]:
                print(f"  [{e.kind.value}] {e.id} — {e.title}")
        else:
            print(f"  result: {r.result}")
    return 0


def _cmd_export(args) -> int:
    m = default_yu_cook_manifest()
    if args.format == "json":
        out = {
            "pi": m.pi, "co_pi": m.co_pi,
            "entries": [_entry_to_dict(e) for e in m.entries],
        }
        json.dump(out, sys.stdout, ensure_ascii=False, indent=2)
        sys.stdout.write("\n")
    else:
        print(f"# Yu-Cook Lab Asset Manifest\n")
        print(f"**PI:** {m.pi}  ")
        print(f"**Co-PI:** {m.co_pi}  ")
        print(f"**Entries:** {len(m)}  \n")
        for kind in EntryKind:
            entries = m.by_kind(kind)
            if not entries:
                continue
            print(f"## {kind.value.upper()} ({len(entries)})\n")
            for e in entries:
                print(f"### {e.title}")
                print(f"- **id:** `{e.id}`")
                print(f"- **status:** {e.status}")
                if e.tags:
                    print(f"- **tags:** {', '.join(e.tags)}")
                if e.n_samples is not None:
                    print(f"- **n:** {e.n_samples}")
                if e.path:
                    print(f"- **path:** `{e.path}`")
                if e.journal_target:
                    print(f"- **target:** {e.journal_target}")
                if e.summary:
                    print(f"\n  {e.summary}")
                print()
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="cancer_discovery_os")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("manifest", help="show lab asset manifest summary")
    p.add_argument("--verbose", "-v", action="store_true")
    p.set_defaults(func=_cmd_manifest)

    p = sub.add_parser("search", help="search lab manifest")
    p.add_argument("term")
    p.add_argument("--verbose", "-v", action="store_true")
    p.set_defaults(func=_cmd_search)

    p = sub.add_parser("ask", help="dispatch a free-text query through the OS")
    p.add_argument("text")
    p.set_defaults(func=_cmd_ask)

    p = sub.add_parser("export", help="export manifest as markdown or json")
    p.add_argument("--format", choices=["markdown", "json"], default="markdown")
    p.set_defaults(func=_cmd_export)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
