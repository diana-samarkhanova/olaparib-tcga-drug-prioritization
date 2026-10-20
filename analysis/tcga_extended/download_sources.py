#!/usr/bin/env python3
"""Retrieve recorded TCGA source snapshots without accepting changed content.

Only entries used in the analysis are retrieved by default. Live source URLs
can change; a SHA-256 mismatch stops the download and preserves any local file.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parent


def fetch(entry):
    if entry.get('retrieval_method') == 'POST retry':
        detail = entry['retry_request']
        req = urllib.request.Request(
            detail['url'], data=json.dumps(detail['body']).encode(),
            headers={'Content-Type': 'application/json', 'Accept': 'application/json'},
        )
        with urllib.request.urlopen(req, timeout=180) as response:
            records = json.loads(response.read())
        gene_id = int(urllib.parse.parse_qs(
            urllib.parse.urlparse(entry['url']).query)['entrezGeneId'][0])
        return json.dumps([x for x in records if x['entrezGeneId'] == gene_id]).encode()
    with urllib.request.urlopen(entry['url'], timeout=180) as response:
        return response.read()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--include-unused', action='store_true',
                        help='Also retrieve provenance-only sources not used in the analysis.')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'source_snapshot_manifest.json').read_text())
    for entry in manifest:
        if not entry.get('used_in_analysis', False) and not args.include_unused:
            continue
        path = ROOT / entry['file']
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256']:
            print('Verified local snapshot:', entry['file'], flush=True)
            continue
        for attempt in range(4):
            try:
                raw = fetch(entry)
                break
            except Exception as error:
                if attempt == 3:
                    raise
                print('Retry:', entry['file'], str(error), flush=True)
                time.sleep(2 ** attempt)
        actual = hashlib.sha256(raw).hexdigest()
        if actual != entry['sha256']:
            raise RuntimeError(
                f"Source snapshot changed for {entry['file']}: expected {entry['sha256']}, "
                f"received {actual}. Local files were preserved. Inspect the source "
                'change before creating a separately documented updated analysis.'
            )
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as output:
                temporary = Path(output.name)
                output.write(raw)
            os.replace(temporary, path)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
        print('Downloaded and verified:', entry['file'], flush=True)


if __name__ == '__main__':
    main()
