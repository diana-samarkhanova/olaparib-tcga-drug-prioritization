#!/usr/bin/env python3
"""Restore local data from Complete_Analysis_EN.zip without replacing source code.

This is an archive replay helper. It restores both inputs and saved results;
restoration itself does not recompute any result. Existing files are retained.
"""
import argparse
import hashlib
import shutil
import stat
import zipfile
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
CODE_SUFFIXES = {'.py', '.pyc', '.r', '.mjs', '.js', '.sh'}


def target_path(name, root):
    parts = PurePosixPath(name).parts
    if not parts or name.startswith('/') or '..' in parts or '\\' in name:
        raise ValueError(f'Unsafe archive path: {name}')
    if parts[0] == 'reproducibility':
        relative = Path(*parts[1:])
        if not relative.parts or relative.parts[0] not in {'analysis', 'upload', 'outputs'}:
            return None
    elif parts[0] == 'tables':
        relative = Path(*parts)
    elif parts[0] == 'figures':
        relative = Path('outputs/TCGA_Drug_Analysis_EN_20261013', *parts)
    else:
        return None
    if relative.suffix.lower() in CODE_SUFFIXES or relative.name == 'requirements.txt':
        return None
    target = (root / relative).resolve()
    if not target.is_relative_to(root):
        raise ValueError(f'Path escapes destination: {name}')
    return target


def restore(archive, root):
    restored = skipped = 0
    with zipfile.ZipFile(archive) as source:
        # Validate every archive entry before creating any file.
        entries = []
        for info in source.infolist():
            if stat.S_ISLNK(info.external_attr >> 16):
                raise ValueError(f'Symlink in archive: {info.filename}')
            target = target_path(info.filename, root)
            if target is not None and not info.is_dir():
                entries.append((info, target))
        for info, target in entries:
            if target.exists():
                skipped += 1
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with source.open(info) as src, target.open('xb') as dst:
                shutil.copyfileobj(src, dst)
            restored += 1
    return restored, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', required=True, type=Path)
    parser.add_argument('--destination', type=Path, default=ROOT)
    args = parser.parse_args()
    root = args.destination.resolve()
    root.mkdir(parents=True, exist_ok=True)
    archive = args.archive.resolve()
    with archive.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    print(f'Archive SHA-256: {digest}')
    restored, skipped = restore(archive, root)
    print(f'Restored {restored} data/result files; retained {skipped} existing files.')
    print('Source code was not replaced. These are restored snapshots, not a new analysis.')
    print('Next: python analysis/tcga_extended/reconstruct_archived_purity.py')


if __name__ == '__main__':
    main()
