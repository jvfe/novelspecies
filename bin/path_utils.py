"""Resolve input file paths across Nextflow work directories."""

from __future__ import annotations

from pathlib import Path


def resolve_input_path(path_str: str, *anchors: Path) -> Path:
    """Return an existing file path, searching relative to anchor directories."""
    raw = Path(path_str)
    if raw.is_file():
        return raw.resolve()

    candidates = [raw.resolve()]
    for anchor in anchors:
        if anchor is None:
            continue
        candidates.append((anchor / raw).resolve())
        if raw.is_absolute():
            continue
        # Also try resolving from project root when paths omit leading directories.
        candidates.append((anchor / raw.name).resolve())

    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if candidate.is_file():
            return candidate

    anchor_list = ", ".join(str(a) for a in anchors if a is not None)
    raise FileNotFoundError(f"Could not resolve input path '{path_str}' (searched anchors: {anchor_list})")
