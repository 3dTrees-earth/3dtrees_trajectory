"""Publish outputs atomically without loading the numerical libraries."""
from __future__ import annotations

import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n", dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(text)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_empty(path: Path) -> None:
    write_text(path, "")


def write_report(path: Path, data: dict) -> None:
    write_text(path, json.dumps(data, indent=1, allow_nan=False) + "\n")
