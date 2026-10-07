"""Rename every file ending in 'updated.csv' to end in '.csv' instead,
in a folder and all its subfolders.

    python strip_updated.py D:/path/to/folder

e.g. 'patient1updated.csv' -> 'patient1.csv'. A file is skipped (and
reported) if the new name already exists, so nothing is overwritten.
"""
import sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
suffix = "_.csv"

for f in root.rglob("*" + suffix):
    if not f.is_file():
        continue
    new = f.with_name(f.name[: -len(suffix)] + ".csv")
    if new.exists():
        print(f"SKIPPED (exists): {f} -> {new.name}")
        continue
    f.rename(new)
    print(f"{f} -> {new.name}")
