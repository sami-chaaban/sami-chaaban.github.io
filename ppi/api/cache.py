from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Optional
import threading
import sys
import time


@dataclass
class CacheEntry:
    value: Any
    created_at: float
    size_bytes: int = 0


class ReportCache:
    def __init__(self, ttl_seconds: int = 86400, max_entries: int = 256,
                 max_bytes: Optional[int] = None, size_of: Optional[Callable[[Any], int]] = None) -> None:
        self.ttl_seconds = ttl_seconds
        self.max_entries = max_entries
        self.max_bytes = max(0, int(max_bytes)) if max_bytes is not None else None
        self._size_of = size_of or self._payload_size
        self._total_bytes = 0
        self._entries: Dict[str, CacheEntry] = {}
        self._lock = threading.Lock()

    def get(self, key: str) -> Optional[Any]:
        with self._lock:
            entry = self._entries.get(key)
            if not entry:
                return None
            if (time.time() - entry.created_at) > self.ttl_seconds:
                self._remove_locked(key)
                return None
            return entry.value

    @staticmethod
    def _payload_size(value: Any) -> int:
        """Account serialized payloads exactly; allow other caches a custom sizer."""
        if isinstance(value, (bytes, bytearray, memoryview)):
            return len(value)
        if isinstance(value, str):
            return len(value.encode('utf-8'))
        size = getattr(value, 'cache_size_bytes', None)
        return int(size) if size is not None else sys.getsizeof(value)

    @property
    def total_bytes(self) -> int:
        with self._lock:
            return self._total_bytes

    def set(self, key: str, value: Any) -> None:
        size = max(0, int(self._size_of(value))) if self.max_bytes is not None else 0
        with self._lock:
            self._remove_locked(key)
            # A single oversized report is returned normally but is not retained.
            if self.max_entries <= 0 or (self.max_bytes is not None and size > self.max_bytes):
                return
            while self._entries and (len(self._entries) >= self.max_entries or
                    (self.max_bytes is not None and self._total_bytes + size > self.max_bytes)):
                self._evict_oldest_locked()
            self._entries[key] = CacheEntry(value=value, created_at=time.time(), size_bytes=size)
            self._total_bytes += size

    def _remove_locked(self, key: str) -> None:
        entry = self._entries.pop(key, None)
        if entry is not None:
            self._total_bytes -= entry.size_bytes

    def _evict_oldest_locked(self) -> None:
        if not self._entries:
            return
        oldest_key = min(self._entries.items(), key=lambda item: item[1].created_at)[0]
        self._remove_locked(oldest_key)
