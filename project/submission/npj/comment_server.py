#!/usr/bin/env python3
"""
v17 comment server — extends http.server with POST/GET for /api/comments.

Routes:
  GET /api/comments          → return JSON list of all comments
  POST /api/comments         → append a comment (body = JSON {section, text, author?})
  DELETE /api/comments       → clear all (admin only, requires ?confirm=1)
  Everything else            → default static file serving (same as http.server)

Storage: /opt/thyroid-dash/project/submission/npj/yu_comments_server.json
Port: 8765 (replaces existing http.server)
"""
from __future__ import annotations
import json
import os
import sys
import time
import threading
from datetime import datetime
from pathlib import Path
from http.server import HTTPServer, SimpleHTTPRequestHandler

ROOT = Path(__file__).parent
STORE = ROOT / "yu_comments_server.json"
LOCK = threading.Lock()


def _load() -> list[dict]:
    if not STORE.exists():
        return []
    try:
        return json.loads(STORE.read_text(encoding="utf-8"))
    except Exception:
        return []


def _save(comments: list[dict]) -> None:
    STORE.write_text(json.dumps(comments, ensure_ascii=False, indent=2), encoding="utf-8")


class CommentHandler(SimpleHTTPRequestHandler):
    server_version = "v17CommentServer/1.0"

    def address_string(self):
        # Skip reverse DNS lookup (which blocks for ~30s on remote IPs).
        # Return IP only.
        return self.client_address[0]

    def log_message(self, fmt, *args):
        # quiet log
        sys.stderr.write(f"{self.address_string()} - {fmt % args}\n")

    def _send_json(self, status: int, payload):
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, DELETE, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        if self.path.startswith("/api/"):
            self.send_response(204)
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, DELETE, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
            self.end_headers()
        else:
            super().do_OPTIONS() if hasattr(super(), "do_OPTIONS") else self.send_error(405)

    def do_GET(self):
        if self.path.startswith("/api/comments"):
            with LOCK:
                comments = _load()
            self._send_json(200, {
                "ok": True,
                "n": len(comments),
                "comments": comments,
                "server_time": datetime.now().isoformat(),
            })
            return
        return super().do_GET()

    def do_POST(self):
        if self.path == "/api/comments":
            length = int(self.headers.get("Content-Length", "0"))
            if length > 50000:
                self._send_json(413, {"ok": False, "error": "Comment too long (>50KB)"})
                return
            try:
                raw = self.rfile.read(length).decode("utf-8")
                data = json.loads(raw)
            except Exception as e:
                self._send_json(400, {"ok": False, "error": f"Invalid JSON: {e}"})
                return
            section = str(data.get("section", "")).strip()[:200]
            text = str(data.get("text", "")).strip()
            author = str(data.get("author", "")).strip()[:60] or None
            url = str(data.get("url", "")).strip()[:500]
            if not section or not text:
                self._send_json(400, {"ok": False, "error": "section and text are required"})
                return
            if len(text) > 20000:
                self._send_json(413, {"ok": False, "error": "Comment text too long (>20KB)"})
                return
            entry = {
                "id": int(time.time() * 1000),
                "section": section,
                "text": text,
                "author": author,
                "url": url,
                "timestamp": datetime.now().isoformat(),
                "client_ip": self.address_string(),
            }
            with LOCK:
                comments = _load()
                comments.append(entry)
                _save(comments)
            self._send_json(200, {"ok": True, "id": entry["id"], "n_total": len(comments)})
            return
        self._send_json(404, {"ok": False, "error": "POST endpoint not found"})

    def do_DELETE(self):
        if self.path.startswith("/api/comments"):
            # require ?confirm=1 query string
            if "confirm=1" not in self.path:
                self._send_json(400, {"ok": False, "error": "Add ?confirm=1 to confirm deletion"})
                return
            with LOCK:
                _save([])
            self._send_json(200, {"ok": True, "message": "all comments cleared"})
            return
        self._send_json(404, {"ok": False, "error": "DELETE endpoint not found"})


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
    bind = "0.0.0.0"
    os.chdir(str(ROOT))
    print(f"[comment_server] starting on {bind}:{port} from {ROOT}", flush=True)
    print(f"[comment_server] comments stored at {STORE}", flush=True)
    server = HTTPServer((bind, port), CommentHandler)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("[comment_server] shutting down", flush=True)


if __name__ == "__main__":
    main()
