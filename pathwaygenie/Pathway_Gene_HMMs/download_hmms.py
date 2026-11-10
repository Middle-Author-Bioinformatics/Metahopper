#!/usr/bin/env python3
"""Download this model panel (or selected genes) using only Python's standard library."""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.request


def build_libraries(models, root):
    """Use accessions as unique NAMEs in bundles; preserve original individual files."""
    folder = root / 'libraries'
    folder.mkdir(parents=True, exist_ok=True)
    groups = {
        'pathway_gene_models.hmm': lambda m: m['group'] == 'function_models',
        'context_ncbi.hmm': lambda m: m['group'] == 'context_only' and m['source_database'] != 'Pfam',
        'context_pfam.hmm': lambda m: m['source_database'] == 'Pfam',
    }
    for filename, select in groups.items():
        text = []
        for model in models:
            if select(model):
                s = (root / model['relative_file']).read_text()
                s = re.sub(r'^NAME\s+.+$', 'NAME  ' + model['accession'], s, count=1, flags=re.M)
                text.append(s.rstrip() + '\n')
        (folder / filename).write_text(''.join(text))


def fetch(model, root, force=False):
    p = root / model['relative_file']
    p.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(3):
        try:
            if p.exists() and not force:
                data = p.read_bytes()
            else:
                request = urllib.request.Request(model['url'], headers={'User-Agent': 'PathwayGeneHMMDownloader/1.0'})
                with urllib.request.urlopen(request, timeout=60) as response:
                    data = response.read()
                if data.startswith(b'\x1f\x8b'):
                    data = gzip.decompress(data)
            if hashlib.sha256(data).hexdigest() != model['sha256']:
                raise ValueError('Checksum mismatch: file is incomplete or the source model has changed. '
                                 'Do not silently substitute a new model version.')
            if not data.startswith(b'HMMER3/'):
                raise ValueError('Response is not a HMMER3 model')
            p.write_bytes(data)
            print('Verified', model['accession'], flush=True)
            return
        except Exception:
            if attempt == 2:
                raise
            time.sleep(attempt + 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--genes', nargs='+', help='Optional case-sensitive target labels, e.g. trpA trpB OAT')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--force', action='store_true', help='Redownload even if local files exist')
    args = parser.parse_args()
    manifest = json.loads((Path(__file__).resolve().parent / 'download_manifest.json').read_text())
    models = manifest['models']
    if args.genes:
        known = {g for m in models for g in m['targets']}
        if set(args.genes) - known:
            parser.error('Unknown target labels: ' + ', '.join(sorted(set(args.genes) - known)))
        models = [m for m in models if set(m['targets']) & set(args.genes)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, min(args.threads, 8))) as pool:
        list(pool.map(lambda m: fetch(m, args.output, args.force), models))
    build_libraries(models, args.output)
    print(f'Complete: {len(models)} verified models. Libraries: {args.output / "libraries"}')


if __name__ == '__main__':
    main()
