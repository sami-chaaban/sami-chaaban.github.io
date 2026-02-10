#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
API_HOST="${API_HOST:-0.0.0.0}"
API_PORT="${API_PORT:-8000}"
WEB_HOST="${WEB_HOST:-0.0.0.0}"
WEB_PORT="${WEB_PORT:-3000}"

cleanup() {
  echo "\nStopping dev servers..."
  if [[ -n "${API_PID:-}" ]]; then
    kill "$API_PID" 2>/dev/null || true
  fi
  if [[ -n "${WEB_PID:-}" ]]; then
    kill "$WEB_PID" 2>/dev/null || true
  fi
}
trap cleanup INT TERM

(
  cd "$ROOT/api"
  echo "Starting API on http://localhost:${API_PORT}"
  python -m uvicorn main:app --reload --host "$API_HOST" --port "$API_PORT"
) &
API_PID=$!

(
  cd "$ROOT/web"
  echo "Starting web on http://localhost:${WEB_PORT}"
  python -m http.server "$WEB_PORT" --bind "$WEB_HOST"
) &
WEB_PID=$!

echo "\nOpen: http://localhost:${WEB_PORT} (API: http://localhost:${API_PORT})"
wait
