from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional
import threading
import time


@dataclass
class CacheEntry:
    value: Any
    created_at: float


class ReportCache:
    def __init__(self, ttl_seconds: int = 86400, max_entries: int = 256) -> None:
        self.ttl_seconds = ttl_seconds
        self.max_entries = max_entries
        self._entries: Dict[str, CacheEntry] = {}
        self._lock = threading.Lock()

    def get(self, key: str) -> Optional[Any]:
        with self._lock:
            entry = self._entries.get(key)
            if not entry:
                return None
            if (time.time() - entry.created_at) > self.ttl_seconds:
                self._entries.pop(key, None)
                return None
            return entry.value

    def set(self, key: str, value: Any) -> None:
        with self._lock:
            if len(self._entries) >= self.max_entries:
                self._evict_oldest_locked()
            self._entries[key] = CacheEntry(value=value, created_at=time.time())

    def _evict_oldest_locked(self) -> None:
        if not self._entries:
            return
        oldest_key = min(self._entries.items(), key=lambda item: item[1].created_at)[0]
        self._entries.pop(oldest_key, None)
