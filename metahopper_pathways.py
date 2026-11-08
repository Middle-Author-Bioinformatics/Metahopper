#!/usr/bin/env python3
"""Run or import PathwayGenie for one MetaHopper rank, preserving verified caches."""
from __future__ import annotations
import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from metahopper_pathway_view import read_results, safe_name, safe_path, sha256

HERE = Path(__file__).resolve().parent
MANAGED = {'-bin_dir', '--bin-dir', '-bin_ext', '--bin-ext', '-out', '--out', '-t', '--threads',
           '-j', '--jobs', '--hmm-dir', '--config', '--force', '--recursive'}


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix='.' + path.name, delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(data, handle, indent=2)
        handle.write('\n')
    temporary.replace(path)


def utc():
    return datetime.now(timezone.utc).isoformat()


def resolve_config(path):
    path = Path(path).expanduser().resolve()
    config = json.loads(path.read_text())
    for key in ('pathwaygenie_script', 'pathwaygenie_hmm_dir', 'pathwaygenie_config'):
        if config.get(key):
            value = Path(os.path.expandvars(config[key])).expanduser()
            config[key] = str((path.parent / value).resolve())
    return config


def runtime_config(config, threads=16, jobs=3):
    if threads < 1 or not 1 <= jobs <= threads:
        raise ValueError('PathwayGenie requires 1 <= jobs <= threads')
    script = Path(config.get('pathwaygenie_script') or HERE / 'pathwaygenie/PathwayGenie.py').expanduser().resolve()
    hmm = Path(config.get('pathwaygenie_hmm_dir') or '').expanduser().resolve()
    rules = Path(config.get('pathwaygenie_config') or HERE / 'pathwaygenie/pathwaygenie_config.json').expanduser().resolve()
    if not script.is_file():
        raise ValueError('PathwayGenie script not found: ' + str(script))
    if not config.get('pathwaygenie_hmm_dir') or not (hmm / 'download_manifest.json').is_file():
        raise ValueError('Set pathwaygenie_hmm_dir to the existing Pathway_Gene_HMMs directory (including download_manifest.json)')
    if not rules.is_file():
        raise ValueError('PathwayGenie pathway configuration not found: ' + str(rules))
    prefix = config.get('tools', {}).get('pathwaygenie')
    if prefix is None:
        env = config.get('pathwaygenie_env', 'pathwaygenie')
        prefix = ['conda', 'run', '--no-capture-output', '-n', env, 'python', str(script)] if env else [sys.executable, str(script)]
    if not isinstance(prefix, list) or not prefix or not all(isinstance(x, str) and x for x in prefix):
        raise ValueError('tools.pathwaygenie must be a complete executable/script prefix array')
    if shutil.which(prefix[0]) is None:
        raise ValueError('PathwayGenie launcher not found: ' + prefix[0])
    extra = config.get('pathwaygenie_args', [])
    if not isinstance(extra, list) or not all(isinstance(x, str) for x in extra):
        raise ValueError('pathwaygenie_args must be an array of strings')
    if any(x.split('=', 1)[0] in MANAGED for x in extra):
        raise ValueError('pathwaygenie_args cannot override managed input, output, configuration, resources or force options')
    return prefix, hmm, rules, extra


def preflight_runtime(config, threads=16, jobs=3):
    prefix, hmm, rules, extra = runtime_config(config, threads, jobs)
    try:
        subprocess.run(prefix + ['--version'], text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, timeout=60, check=True)
        if len(prefix) >= 2 and prefix[-1].endswith('.py') and Path(prefix[-2]).name.startswith('python'):
            tools = {'prodigal': 'prodigal', 'hmmsearch': 'hmmsearch'}
            for i, token in enumerate(extra):
                name = token.split('=', 1)[0]
                if name in ('--prodigal', '--hmmsearch'):
                    tools[name[2:]] = token.split('=', 1)[1] if '=' in token else extra[i + 1]
            probe = ('import shutil,sys; names=' + repr(tools) + '; missing=[k+": "+v for k,v in names.items() if not shutil.which(v)]; '
                     'print("Missing PathwayGenie tools: "+", ".join(missing)) if missing else None; sys.exit(bool(missing))')
            subprocess.run(prefix[:-1] + ['-c', probe], text=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=60, check=True)
    except (subprocess.SubprocessError, OSError, IndexError) as exc:
        detail = getattr(exc, 'stdout', '') or str(exc)
        if isinstance(detail, bytes):
            detail = detail.decode(errors='replace')
        raise ValueError('PathwayGenie environment/dependency check failed: ' + detail[-2000:]) from exc
    return prefix, hmm, rules, extra


def target_directory(root, rank='genus', stage='final'):
    root = Path(root).resolve()
    rank = safe_name(rank)
    base = root / 'final' if stage == 'final' else root
    if not (base / 'bins' / rank).is_dir():
        raise ValueError('No ' + stage + ' bins at rank ' + rank + ': ' + str(base / 'bins' / rank))
    for path in (base / 'bins' / rank).glob('*.fasta'):
        safe_path(path, root)
    return base / 'pathways' / rank


@contextmanager
def result_lock(directory):
    directory.parent.mkdir(parents=True, exist_ok=True)
    path = directory.parent / ('.' + directory.name + '.metahopper-pathways.lock')
    with path.open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('Another pathway run/import is using ' + str(directory)) from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def run_pathways(root, config=None, rank='genus', stage='final', threads=16, jobs=3, force=False, log=print):
    root = Path(root).resolve()
    config = config or {}
    output = target_directory(root, rank, stage)
    prefix, hmm, rules, extra = preflight_runtime(config, threads, jobs)
    base = root / 'final' if stage == 'final' else root
    command = prefix + ['--bin-dir', str(base / 'bins' / rank), '--bin-ext', 'fasta', '--out', str(output),
                        '--hmm-dir', str(hmm), '--config', str(rules), '-t', str(threads), '-j', str(jobs)] + extra
    if force:
        command.append('--force')
    logfile = root / 'logs' / ('pathwaygenie.' + rank + '.log')
    statusfile = root / ('pathwaygenie.' + rank + '.status.json')
    logfile.parent.mkdir(parents=True, exist_ok=True)
    record = dict(schema=1, module='pathwaygenie', rank=rank, stage=stage, status='running',
                  started_utc=utc(), output=output.relative_to(root).as_posix(), command=command,
                  log=logfile.relative_to(root).as_posix())
    started = time.monotonic()
    with result_lock(output):
        atomic_json(statusfile, record)
        log('PathwayGenie: ' + rank + ' bins; verified Prodigal/HMMER caches will be reused. Log: ' + str(logfile))
        try:
            with logfile.open('a') as handle:
                handle.write('\nSTART ' + record['started_utc'] + '\nCOMMAND ' + json.dumps(command) + '\n')
                handle.flush()
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
            if result.returncode not in (0, 1):
                raise RuntimeError('PathwayGenie exited ' + str(result.returncode) + '; see ' + str(logfile))
            data = read_results(output, root, stage, rank, check_wrapper=False)
            failed = [b['bin_id'] for b in data['bins'] if b['status'] == 'failed']
            if result.returncode and not failed:
                raise RuntimeError('PathwayGenie returned an error without finished failed-bin records; see ' + str(logfile))
            record.update(status='partial' if failed else 'completed', failed_bins=failed,
                          exit_code=result.returncode, finished_utc=utc(), elapsed_seconds=round(time.monotonic() - started, 2))
            atomic_json(statusfile, record)
            log('PathwayGenie ' + record['status'] + ': ' + str(output))
            return dict(status=record['status'], failed_bins=failed, output=output,
                        manifest=output / 'run_manifest.json', warnings=data['warnings'])
        except BaseException as exc:
            record.update(status='failed', error=str(exc), finished_utc=utc(), elapsed_seconds=round(time.monotonic() - started, 2))
            atomic_json(statusfile, record)
            raise


def import_results(root, source, rank='genus', stage='final', force=False):
    root, source = Path(root).resolve(), Path(source).expanduser().resolve()
    output = target_directory(root, rank, stage)
    if not (source / 'pathway_summary.long.csv').is_file():
        raise ValueError('Select the extracted PathwayGenie output directory containing pathway_summary.long.csv')
    if (source / '.pathwaygenie.lock').exists():
        raise ValueError('Source PathwayGenie output is locked; finish its run before importing')
    if source == output:
        return read_results(output, root, stage, rank)
    if source.is_relative_to(output) or output.is_relative_to(source):
        raise ValueError('Import source and destination cannot contain one another')
    with result_lock(output):
        if output.exists() and not force:
            raise ValueError(str(output) + ' already exists; use --force-overwrite to import over it (a backup is kept)')
        original_manifest = sha256(source / 'run_manifest.json')
        temporary = Path(tempfile.mkdtemp(prefix='.pathway-import-', dir=output.parent))
        backup = None
        try:
            for path in source.rglob('*'):
                relative = path.relative_to(source)
                if any(part.startswith('.') for part in relative.parts):
                    continue
                if path.is_symlink():
                    raise ValueError('Import source contains a symlink: ' + str(path))
                if path.is_file():
                    destination = temporary / relative
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(path, destination)
            if (source / '.pathwaygenie.lock').exists() or sha256(source / 'run_manifest.json') != original_manifest:
                raise ValueError('Source results changed during import; retry after the run finishes')
            data = read_results(temporary, root, stage, rank, check_wrapper=False)
            if not any(b['source_state'] == 'current' for b in data['bins']):
                raise ValueError('No PathwayGenie input checksums match this sample/rank; results were not imported')
            if output.exists():
                backup = output.with_name(output.name + '.before-import-' + str(time.time_ns()))
                output.replace(backup)
            try:
                temporary.replace(output)
            except BaseException:
                if backup is not None and not output.exists():
                    backup.replace(output)
                raise
            atomic_json(root / ('pathwaygenie.' + rank + '.status.json'), dict(schema=1,module='pathwaygenie',
                        status='imported',rank=rank,stage=stage,finished_utc=utc(),output=output.relative_to(root).as_posix()))
        finally:
            if temporary.exists():
                shutil.rmtree(temporary)
    return read_results(output, root, stage, rank)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['run', 'import'])
    parser.add_argument('-i', '--input', type=Path, required=True, help='MetaHopper sample run directory')
    parser.add_argument('-c', '--config', type=Path, help='MetaHopper config.json (optional for standalone runs)')
    parser.add_argument('--rank', default='genus')
    parser.add_argument('--bin-set', choices=['final', 'preliminary'], default='final')
    parser.add_argument('--results-dir', type=Path, help='Existing extracted PathwayGenie output folder to import')
    parser.add_argument('--pathwaygenie-dir', type=Path, help='Existing installation with PathwayGenie.py, rules and Pathway_Gene_HMMs')
    parser.add_argument('--hmm-dir', type=Path)
    parser.add_argument('--pathway-config', type=Path)
    parser.add_argument('--conda-env', help='Environment containing PathwayGenie dependencies; default pathwaygenie')
    parser.add_argument('-t', '--threads', type=int, default=16)
    parser.add_argument('-j', '--jobs', type=int, default=3)
    parser.add_argument('--force-overwrite', action='store_true', help='Run: rerun external stages. Import: retain a backup and replace existing results.')
    args = parser.parse_args(argv)
    try:
        if args.action == 'import':
            if args.results_dir is None:
                raise ValueError('import requires --results-dir')
            data = import_results(args.input, args.results_dir, args.rank, args.bin_set, args.force_overwrite)
            print('Imported ' + str(len(data['bins'])) + ' bins into ' + data['directory'])
            return 1 if any(b['status'] == 'failed' for b in data['bins']) else 0
        config = resolve_config(args.config) if args.config else {}
        if args.pathwaygenie_dir:
            directory = args.pathwaygenie_dir.expanduser().resolve()
            config.update(pathwaygenie_script=str(directory / 'PathwayGenie.py'),
                          pathwaygenie_hmm_dir=str(directory / 'Pathway_Gene_HMMs'),
                          pathwaygenie_config=str(directory / 'pathwaygenie_config.json'))
        if args.hmm_dir:
            config['pathwaygenie_hmm_dir'] = str(args.hmm_dir.expanduser().resolve())
        if args.pathway_config:
            config['pathwaygenie_config'] = str(args.pathway_config.expanduser().resolve())
        if args.conda_env is not None:
            config['pathwaygenie_env'] = args.conda_env
        result = run_pathways(args.input, config, args.rank, args.bin_set, args.threads, args.jobs, args.force_overwrite)
        return 1 if result['status'] != 'completed' else 0
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError) as exc:
        parser.exit(1, 'ERROR: ' + str(exc) + '\n')


if __name__ == '__main__':
    raise SystemExit(main())
