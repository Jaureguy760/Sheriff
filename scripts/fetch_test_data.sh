#!/usr/bin/env bash
set -euo pipefail

# Fetch the chr19 test tarball (used for quick perf/validation runs).
# Default URL points to the GitHub release asset; override with TEST_DATA_URL if needed.

DEFAULT_URL="https://github.com/BradBalderson/Sheriff/releases/download/test-data-v1/test_data_chr19.tar.gz"
URL="${TEST_DATA_URL:-$DEFAULT_URL}"
OUT="test_data_chr19.tar.gz"
SHA="49e29735a7667c395f0d85c789134fcdf4a21cdb14929b780ddf49237f0fa91d"

echo "[fetch_test_data] target: $OUT"
echo "[fetch_test_data] url:    $URL"

if [[ -f "$OUT" ]]; then
  echo "[fetch_test_data] file exists, verifying checksum..."
  if echo "$SHA  $OUT" | sha256sum --check --status; then
    echo "[fetch_test_data] checksum OK, nothing to do."
    exit 0
  else
    echo "[fetch_test_data] checksum mismatch, re-downloading..."
    rm -f "$OUT"
  fi
fi

curl -L "$URL" -o "$OUT"
echo "[fetch_test_data] downloaded, verifying checksum..."
if echo "$SHA  $OUT" | sha256sum --check --status; then
  echo "[fetch_test_data] success."
else
  echo "[fetch_test_data] checksum FAILED; remove $OUT and retry (or set TEST_DATA_URL)." >&2
  exit 1
fi
