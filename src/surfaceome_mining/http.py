"""Shared HTTP client for public surfaceome APIs.

Protocol:
- identify the client with a User-Agent
- accept gzip
- retry 429 using Retry-After (default 60s) and retry transient 5xx
- pace per-gene calls; prefer bulk catalog/TSV endpoints for genome-wide work
- treat X-Surfaceome-Degraded as a stale response and retry once
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import httpx

DEFAULT_USER_AGENT = (
    "surfaceomeMining/0.1 (+https://github.com/mightenyip/surfaceomeMining)"
)


class SourceHttpError(RuntimeError):
    def __init__(self, message: str, *, status: int | None = None, body: Any = None):
        super().__init__(message)
        self.status = status
        self.body = body


class PacedClient:
    def __init__(
        self,
        *,
        user_agent: str = DEFAULT_USER_AGENT,
        min_interval_s: float = 0.4,
        timeout_s: float = 60.0,
        cache_dir: Path | None = None,
        max_retries: int = 4,
    ) -> None:
        self._client = httpx.Client(
            headers={
                "User-Agent": user_agent,
                "Accept": "application/json, text/markdown, text/tab-separated-values, */*",
                "Accept-Encoding": "gzip, deflate",
            },
            timeout=timeout_s,
            follow_redirects=True,
        )
        self.min_interval_s = min_interval_s
        self.cache_dir = cache_dir
        self.max_retries = max_retries
        self._last_request_at = 0.0

    def close(self) -> None:
        self._client.close()

    def __enter__(self) -> "PacedClient":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def get_json(self, url: str, *, params: dict[str, Any] | None = None, cache_key: str | None = None) -> Any:
        payload, _ = self.get_bytes(url, params=params, cache_key=cache_key)
        if not payload:
            return None
        return json.loads(payload.decode("utf-8"))

    def get_text(self, url: str, *, params: dict[str, Any] | None = None) -> str:
        payload, _ = self.get_bytes(url, params=params)
        return payload.decode("utf-8")

    def get_bytes(
        self,
        url: str,
        *,
        params: dict[str, Any] | None = None,
        cache_key: str | None = None,
    ) -> tuple[bytes, httpx.Headers]:
        cached = self._read_cache(cache_key) if cache_key else None
        if cached is not None:
            return cached, httpx.Headers()

        last_error: Exception | None = None
        for attempt in range(self.max_retries):
            self._pace()
            try:
                response = self._client.get(url, params=params)
            except httpx.HTTPError as exc:
                last_error = exc
                time.sleep(min(2**attempt, 8))
                continue

            if response.status_code == 429:
                wait_s = _retry_after_seconds(response, default=60.0)
                time.sleep(wait_s)
                continue
            if response.status_code >= 500:
                time.sleep(min(2**attempt, 8))
                continue
            if response.status_code == 404:
                raise SourceHttpError(
                    f"not found: {url}",
                    status=404,
                    body=_maybe_json(response),
                )
            if response.status_code >= 400:
                raise SourceHttpError(
                    f"HTTP {response.status_code} for {url}",
                    status=response.status_code,
                    body=_maybe_json(response),
                )
            if response.headers.get("X-Surfaceome-Degraded") and attempt < self.max_retries - 1:
                time.sleep(1.0)
                continue

            content = response.content
            if cache_key:
                self._write_cache(cache_key, content)
            return content, response.headers

        raise SourceHttpError(f"failed after retries: {url}") from last_error

    def _pace(self) -> None:
        elapsed = time.monotonic() - self._last_request_at
        if elapsed < self.min_interval_s:
            time.sleep(self.min_interval_s - elapsed)
        self._last_request_at = time.monotonic()

    def _cache_path(self, cache_key: str) -> Path | None:
        if not self.cache_dir:
            return None
        return self.cache_dir / cache_key

    def _read_cache(self, cache_key: str) -> bytes | None:
        path = self._cache_path(cache_key)
        if path is None or not path.is_file():
            return None
        return path.read_bytes()

    def _write_cache(self, cache_key: str, content: bytes) -> None:
        path = self._cache_path(cache_key)
        if path is None:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)


def _retry_after_seconds(response: httpx.Response, default: float) -> float:
    raw = response.headers.get("Retry-After")
    if not raw:
        return default
    try:
        return max(float(raw), 1.0)
    except ValueError:
        return default


def _maybe_json(response: httpx.Response) -> Any:
    try:
        return response.json()
    except ValueError:
        return response.text[:500]
