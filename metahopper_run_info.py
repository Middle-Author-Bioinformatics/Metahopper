"""Downloadable run settings and complete saved execution logs for one sample."""
import hashlib
import json
from pathlib import Path
import tempfile


def run_files(root, stage):
    root = Path(root).resolve()
    output = root / 'report_support'
    if not output.resolve().is_relative_to(root):
        raise ValueError('Report support directory points outside this run')
    output.mkdir(exist_ok=True)

    def checked(path):
        if not path.resolve(strict=True).is_relative_to(root):
            raise ValueError('Run information or log points outside this run: ' + str(path))
        return path

    records = set()
    for name in ('run.json', 'run_manifest.json', 'pipeline_status.json', 'postprocess_status.json'):
        if (root / name).is_file():
            records.add(root / name)
    records.update(root.glob('pathwaygenie.*.status.json'))
    for base in (root, root / 'final'):
        for pattern in ('pathways/*/run_manifest.json', 'pathways/*/config_used.json',
                        'annotation/*/annotation_manifest.json'):
            records.update(base.glob(pattern))
    settings = {}
    for path in sorted(records):
        checked(path)
        relative = path.relative_to(root).as_posix()
        content = path.read_bytes()
        try:
            value = json.loads(content.decode('utf-8-sig'))
        except (ValueError, UnicodeError) as exc:
            value = dict(read_error=str(exc))
        settings[relative] = dict(sha256=hashlib.sha256(content).hexdigest(), data=value)

    logs = set()
    for pattern in ('pipeline.log*', 'postprocess.log*', 'pathwaygenie.*.log'):
        logs.update(path for path in root.glob(pattern) if path.is_file())
    for directory in ('logs', 'reassembly', 'annotation', 'final/annotation', 'pathways', 'final/pathways'):
        base = root / directory
        if base.is_dir():
            logs.update(path for path in base.rglob('*.log') if path.is_file())
    if (root / 'logs/commands.jsonl').is_file():
        logs.add(root / 'logs/commands.jsonl')
    catalog = []
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=output, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(('MetaHopper saved execution logs\nSource run: ' + root.name
                          + '\nComplete contents of the saved logs listed below; no truncation.\n').encode())
            for path in sorted(logs):
                checked(path)
                before = path.stat()
                relative = path.relative_to(root).as_posix()
                stream.write(('\n\n========== ' + relative + ' ==========\n').encode())
                digest = hashlib.sha256()
                with path.open('rb') as source:
                    for block in iter(lambda: source.read(1024 * 1024), b''):
                        digest.update(block)
                        stream.write(block)
                after = path.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise ValueError('Execution log changed while preparing the report: ' + relative)
                catalog.append(dict(path=relative, size=after.st_size, sha256=digest.hexdigest()))
            if not catalog:
                stream.write(b'\nNo saved execution logs were found for this run.\n')
        temporary.replace(output / 'full_run.log')
        info = dict(schema='metahopper-run-info-1', source_run=root.name, bin_set=stage,
                    settings=settings, logs=catalog, log_file='report_support/full_run.log',
                    note='Saved execution records at report generation; absent stages or logs are not inferred.')
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=output, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(info, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
        temporary.replace(output / 'run_info.json')
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return dict(info='report_support/run_info.json', log='report_support/full_run.log' if catalog else None)
