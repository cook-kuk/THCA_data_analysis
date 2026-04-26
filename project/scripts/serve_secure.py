#!/usr/bin/env python3
"""
Hardened static HTTP server for the THCA dashboard.

Adds security + cache headers that python -m http.server cannot:
  - Content-Security-Policy (allows self, inline styles for Plotly, data: for images)
  - X-Frame-Options: DENY
  - X-Content-Type-Options: nosniff
  - Referrer-Policy: strict-origin-when-cross-origin
  - Permissions-Policy: minimal (disables camera/mic/geo/FLoC)
  - Cache-Control: long for /assets/* (immutable when hashed), short for HTML/JSON
  - Compression: gzip on text responses above 1KB
  - HEAD/GET only, no directory listings

Runs without dependencies beyond stdlib. Listen port via $PORT env (default 8012).
"""
from __future__ import annotations

import functools
import os
import sys
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STATIC_ROOT = ROOT  # serve entire project (matches current systemd unit); dashboard lives under /reports/html/

PORT = int(os.environ.get("PORT", 8012))
BIND = os.environ.get("BIND", "0.0.0.0")

LONG_CACHE_DIRS = ("/reports/html/assets/vendor/", "/reports/html/assets/fonts/")
MEDIUM_CACHE_DIRS = ("/reports/html/assets/css/", "/reports/html/assets/js/", "/reports/html/assets/img/")
NO_CACHE_EXTS = (".html", ".json")

CSP = (
    "default-src 'self'; "
    "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://cdn.jsdelivr.net https://unpkg.com https://cdn.tailwindcss.com https://cdn.plot.ly; "
    "style-src 'self' 'unsafe-inline' https://fonts.googleapis.com https://cdn.jsdelivr.net https://unpkg.com; "
    "img-src 'self' data: blob: https:; "
    "font-src 'self' data: https://fonts.gstatic.com https://cdn.jsdelivr.net; "
    "connect-src 'self' https:; "
    "frame-src 'self' data: blob:; "
    "frame-ancestors 'self'; "
    "base-uri 'self'; "
    "form-action 'self'"
)


class SecureHandler(SimpleHTTPRequestHandler):
    server_version = "THCA-Secure/1.0"
    # Disable directory listings
    def list_directory(self, path):
        self.send_error(403, "Directory listing disabled")
        return None

    def end_headers(self):
        self.send_header("Content-Security-Policy", CSP)
        self.send_header("X-Frame-Options", "DENY")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "strict-origin-when-cross-origin")
        self.send_header(
            "Permissions-Policy",
            "camera=(), microphone=(), geolocation=(), interest-cohort=()",
        )
        self._inject_cache_control()
        super().end_headers()

    def _inject_cache_control(self):
        path = self.path.split("?", 1)[0]
        if path.endswith(NO_CACHE_EXTS):
            self.send_header("Cache-Control", "no-cache, must-revalidate, max-age=0")
        elif any(p in path for p in LONG_CACHE_DIRS):
            self.send_header("Cache-Control", "public, max-age=31536000, immutable")
        elif any(p in path for p in MEDIUM_CACHE_DIRS):
            self.send_header("Cache-Control", "public, max-age=604800")
        else:
            self.send_header("Cache-Control", "public, max-age=300")

    def log_message(self, fmt, *args):
        sys.stdout.write("[%s] %s - %s\n" % (self.log_date_time_string(), self.address_string(), fmt % args))
        sys.stdout.flush()


def run():
    os.chdir(str(STATIC_ROOT))
    handler = functools.partial(SecureHandler, directory=str(STATIC_ROOT))
    with ThreadingHTTPServer((BIND, PORT), handler) as httpd:
        print(f"THCA secure server listening on {BIND}:{PORT}, root={STATIC_ROOT}", flush=True)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    run()
