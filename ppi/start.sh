#!/usr/bin/env bash
set -euo pipefail

PPI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PPI_ROOT"
# Docker has the reviewed patch baked in. Local startup uses a project-owned
# package overlay and does not modify the caller's shared Python environment.
if [[ "${ROAMI_ARPEGGIO_PATCHED:-0}" != "1" ]]; then
  "${CHAPI_PYTHON:-python}" scripts/patch_arpeggio.py --overlay "$PPI_ROOT/.arpeggio/lib"
  export PYTHONPATH="$PPI_ROOT/.arpeggio/lib:$PPI_ROOT${PYTHONPATH:+:$PYTHONPATH}"
fi

if [[ "${RENDER:-}" == "true" ]]; then
  export CHAPI_LOW_MEMORY_MODE="${CHAPI_LOW_MEMORY_MODE:-1}"
  export CHAPI_PERSISTENT_WORKER="${CHAPI_PERSISTENT_WORKER:-0}"
  export CHAPI_WORKER_WARMUP_AT_STARTUP="${CHAPI_WORKER_WARMUP_AT_STARTUP:-0}"
  export CHAPI_SPLIT_CHAIN_BATCH_LIMIT="${CHAPI_SPLIT_CHAIN_BATCH_LIMIT:-1}"
  export CHAPI_MESH_CACHE_MAX_ENTRIES="${CHAPI_MESH_CACHE_MAX_ENTRIES:-1}"
  export CHAPI_MESH_CACHE_MAX_BYTES="${CHAPI_MESH_CACHE_MAX_BYTES:-0}"
  export STRUCTURE_TEXT_CACHE_MAX_ENTRIES="${STRUCTURE_TEXT_CACHE_MAX_ENTRIES:-1}"
fi

exec "${CHAPI_PYTHON:-python}" -m uvicorn api.main:app --host "${API_HOST:-0.0.0.0}" --port "${PORT:-8000}"
