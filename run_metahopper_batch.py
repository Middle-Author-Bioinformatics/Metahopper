#!/usr/bin/env python3
"""Run paired Aphrophora samples sequentially, then publish each completed report.

Run with the Python in metahopper-lean. No additional Python packages required.
Default: discover *_R1/R2.fq.gz or .fastq.gz, excluding APHCOR_AUC2351.
Use --reuse-audit to inspect saved assemblies/searches without running analysis.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from dataclasses import dataclass
from datetime import datetime
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace

BUCKET = 'amplify-d2axfpg5z1hc8g-ma-metahopperreportsbucket7-iiwufem9jwha'
RANKS = ['phylum', 'class', 'order', 'family', 'genus', 'species']
SUMMARY_FIELDS = ['sample', 'output_directory', 'status', 'stage', 'started', 'finished',
                  'detail', 'pipeline_log', 'publish_log']


@dataclass(frozen=True)
class Sample:
    name: str
    read1: Path
    read2: Path
    output: Path


def now():
    return datetime.now().astimezone().isoformat(timespec='seconds')


def saved_reads_match(folder, read1, read2):
    manifest = json.loads((folder / 'run.json').read_text())
    arguments = manifest.get('arguments', {})
    def matches(value, expected):
        if not value:
            return False
        path = Path(os.path.expandvars(str(value))).expanduser()
        possibilities = [path] if path.is_absolute() else [folder.parent / path, folder / path]
        return any(candidate.resolve() == expected.resolve() for candidate in possibilities)
    return matches(arguments.get('read1'), read1) and matches(arguments.get('read2'), read2)


def choose_output(source, name, read1, read2):
    preferred = source / ('APHCOR_AUC2351-2' if name == 'APHCOR_AUC2351' else name)
    candidates = [preferred]
    if name.endswith('_combined'):
        candidates.append(source / name[:-len('_combined')])
    for folder in candidates:
        if (folder / 'run.json').is_file():
            if not saved_reads_match(folder, read1, read2):
                raise ValueError('Saved read pair does not match the selected inputs in ' + str(folder / 'run.json')
                                 + '. Review these inputs before reusing this run.')
            return folder
        if folder.exists() and any(folder.iterdir()):
            raise ValueError('Existing output has no usable run.json: ' + str(folder))
    # Also recognize a renamed run when its saved input pair identifies it.
    matches = []
    for manifest in sorted(source.glob('*/run.json')):
        if manifest.parent in candidates:
            continue
        try:
            if saved_reads_match(manifest.parent, read1, read2):
                matches.append(manifest.parent)
        except (OSError, ValueError, TypeError, AttributeError):
            continue
    if len(matches) > 1:
        raise ValueError('Several saved runs use these inputs: ' + ', '.join(map(str, matches)))
    return matches[0] if matches else preferred


def discover_samples(source, selected=None, excluded=()):
    pairs = {}
    for extension in ('.fq.gz', '.fastq.gz'):
        for mate in (1, 2):
            suffix = '_R' + str(mate) + extension
            for path in source.glob('*' + suffix):
                name = path.name[:-len(suffix)]
                if (selected and name not in selected) or name in excluded:
                    continue
                if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name):
                    raise ValueError('Unsupported sample name: ' + name)
                pair = pairs.setdefault(name, {})
                if mate in pair:
                    raise ValueError('Multiple files for the same read mate: ' + name)
                pair[mate] = path
    if selected and set(selected) - set(pairs):
        raise ValueError('No read pairs for: ' + ', '.join(sorted(set(selected) - set(pairs))))
    samples = []
    for name, pair in sorted(pairs.items()):
        if set(pair) != {1, 2}:
            raise ValueError('Missing R1 or R2 mate for ' + name)
        if any(not p.is_file() or p.stat().st_size == 0 for p in pair.values()):
            raise ValueError('Missing or empty FASTQ: ' + name)
        output = choose_output(source, name, pair[1], pair[2])
        samples.append(Sample(name, pair[1].resolve(), pair[2].resolve(), output))
    if not samples:
        raise ValueError('No selected paired FASTQs found in ' + str(source))
    return samples


def pipeline_command(args, sample):
    command = [sys.executable, str(args.pipeline_dir / 'MetaHopper.py')]
    if (sample.output / 'run.json').is_file():
        command += ['--resume', str(sample.output)]
    else:
        command += ['-1', str(sample.read1), '-2', str(sample.read2), '-o', str(sample.output),
                    '--rescue-eukarya', '--skip-phyloflash']
    command += [
               '-c', str(args.config), '-t', str(args.threads), '--memory', str(args.memory),
               '--checkm',
               '--annotate', 'bakta', '--annotation-rank', 'genus',
               '--pathways', 'pathwaygenie', '--pathway-rank', 'genus',
               '--pathway-jobs', str(args.pathway_jobs)]
    if args.force_overwrite:
        command.append('--force-overwrite')
    return command


def load_pipeline(args):
    sys.path.insert(0, str(args.pipeline_dir))
    spec = importlib.util.spec_from_file_location('_metahopper_batch_preflight', args.pipeline_dir / 'MetaHopper.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, 'parse_args') or not hasattr(module, 'audit_reuse'):
        raise ValueError('Install the updated MetaHopper.py and metahopper_state.py supplied with this resume fix.')
    return module


def load_config(args):
    config = json.loads(args.config.read_text())
    if not isinstance(config, dict):
        raise ValueError('MetaHopper config must be a JSON object: ' + str(args.config))
    for key in ('diamond_db', 'taxdump', 'phyloflash_db', 'bakta_db', 'annotation_metadata',
                'pathwaygenie_script', 'pathwaygenie_hmm_dir', 'pathwaygenie_config'):
        if config.get(key):
            value = Path(os.path.expandvars(config[key])).expanduser()
            config[key] = str((args.config.parent / value).resolve())
    return config


def sample_options(module, args, sample):
    options = module.parse_args(pipeline_command(args, sample)[2:])
    for key in ('read1', 'read2', 'assembly', 'assembly_graph', 'references'):
        if getattr(options, key, None):
            setattr(options, key, str(Path(getattr(options, key)).expanduser().resolve()))
    return options


def audit_samples(args, samples, module=None, config=None, enforce=False):
    module = module or load_pipeline(args)
    config = config if config is not None else load_config(args)
    blockers = []
    for sample in samples:
        print('\n' + sample.name + ' -> ' + str(sample.output), flush=True)
        if not (sample.output / 'run.json').is_file():
            print('  NEW: no saved run; the pipeline will start this sample.', flush=True)
            continue
        options = sample_options(module, args, sample)
        for row in module.audit_reuse(options, config, sample.output):
            print('  ' + row['stage'] + ': ' + row['status'] + ' — ' + row['reason'], flush=True)
            if row['stage'] in ('assembly', 'preliminary_diamond') and row['retained'] and row['status'] != 'REUSE':
                blockers.append(sample.name + ': ' + row['stage'])
    if blockers and enforce and not (args.allow_recompute or args.force_overwrite):
        raise ValueError('Retained assemblies/searches could not be verified for reuse: ' + '; '.join(blockers)
                         + '. No analysis was started. Review --reuse-audit; --allow-recompute explicitly permits rebuilding invalid stages.')
    return blockers


def publish_command(args, sample):
    command = [sys.executable, str(args.publisher), '--source', str(args.source),
               '--sample', sample.output.name, '--pipeline-dir', str(args.pipeline_dir),
               '--bucket', args.bucket, '--prefix', args.prefix,
               '--graph-max-nodes', str(args.graph_max_nodes),
               '--megahit-max-links', str(args.megahit_max_links), '--pathways-rank', 'genus']
    if args.aws_profile:
        command += ['--profile', args.aws_profile]
    return command


def validation_options(args, sample):
    """Options used by MetaHopper.validate, matching pipeline_command above."""
    return SimpleNamespace(read1=str(sample.read1), read2=str(sample.read2), assembly=None,
                           assembly_graph=None, ranks=RANKS.copy(), skip_reassembly=False,
                           reassemble_rank='genus', skip_phyloflash=True, prefilter='tiara',
                           checkm=True, inspect='endosymbionts', skip_report=False,
                           threads=args.threads, memory=args.memory, annotate='bakta',
                           annotation_rank='genus', annotation_timeout_h=6,
                           pathways='pathwaygenie', pathway_rank='genus',
                           pathway_jobs=args.pathway_jobs)


def preflight(args, sample):
    for path in (args.pipeline_dir / 'MetaHopper.py', args.publisher, args.config):
        if not path.is_file():
            raise ValueError('Required file not found: ' + str(path))
    checks = [(args.pipeline_dir / 'MetaHopper.py', ('--pathways', '--annotation-rank')),
              (args.publisher, ('--pathways-rank', '--graph-max-nodes', '--megahit-max-links'))]
    for script, flags in checks:
        result = subprocess.run([sys.executable, str(script), '--help'], text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60, check=True)
        if not all(flag in result.stdout for flag in flags):
            raise ValueError('Install the latest complete Pathways update: ' + str(script))
    if not args.publish_only:
        config = load_config(args)
        module = load_pipeline(args)
        module.validate(sample_options(module, args, sample), config)
        print('MetaHopper, Bakta database, PathwayGenie panel and environment checks passed.', flush=True)
    if shutil.which('aws') is None:
        raise ValueError('AWS CLI is missing from this environment/PATH.')
    command = ['aws'] + (['--profile', args.aws_profile] if args.aws_profile else [])
    subprocess.run(command + ['sts', 'get-caller-identity', '--output', 'json'],
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=60, check=True)
    print('AWS credentials are available. S3 upload permissions are checked when publishing.', flush=True)


def verify_completed(output):
    """A zero exit status alone does not prove that optional analyses succeeded."""
    run = json.loads((output / 'run.json').read_text())
    state = json.loads((output / 'pipeline_status.json').read_text())
    if run.get('status') != 'completed' or state.get('current_step'):
        raise ValueError('Sample pipeline has not completed: ' + str(output))
    steps = state.get('steps', {})
    for name in ('checkm', 'annotation', 'pathways', 'report'):
        record = steps.get(name, {})
        if record.get('status') in ('completed', 'reused'):
            continue
        if name == 'checkm' and record.get('status') == 'skipped' and record.get('message') == 'No named final bins to assess':
            continue
        raise ValueError(name + ' is ' + str(record.get('status', 'missing')) + ': ' + record.get('message', ''))
    report = output / 'metahopper_report.html'
    if not report.is_file() or not report.stat().st_size:
        raise ValueError('Report is missing or empty: ' + str(report))
    summary = output / 'final/annotation/genus/summary.tsv'
    with summary.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or any(row.get('status') not in ('completed', 'skipped') for row in rows):
        raise ValueError('Bakta has failed/stale/missing bin annotations: ' + str(summary))
    pathway = json.loads((output / 'pathwaygenie.genus.status.json').read_text())
    if pathway.get('status') != 'completed':
        raise ValueError('PathwayGenie is not fully completed: ' + str(pathway.get('status')))
    if not (output / 'final/pathways/genus/run_manifest.json').is_file():
        raise ValueError('PathwayGenie results manifest is missing.')


def save_summary(path, rows):
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, newline='', delete=False) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def run_logged(command, logfile):
    print('  ' + shlex.join(command), flush=True)
    print('  Log: ' + str(logfile), flush=True)
    with logfile.open('a') as handle:
        handle.write('\n' + now() + '\n' + shlex.join(command) + '\n')
        handle.flush()
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError('Command exited ' + str(result.returncode) + '; see ' + str(logfile))


@contextmanager
def batch_lock(source):
    with (source / '.metahopper_batch.lock').open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('Another batch runner is active in ' + str(source)) from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def run_batch(args, samples):
    preflight(args, samples[0])
    if not args.publish_only:
        audit_samples(args, samples, enforce=True)
    if args.preflight_only:
        return 0
    base = args.source / 'batch_logs'
    base.mkdir(exist_ok=True)
    logs = Path(tempfile.mkdtemp(prefix=datetime.now().strftime('%Y%m%d_%H%M%S_'), dir=base))
    summary = logs / 'summary.tsv'
    print('Batch summary: ' + str(summary), flush=True)
    rows = []
    for index, sample in enumerate(samples, 1):
        row = dict.fromkeys(SUMMARY_FIELDS, '')
        row.update(sample=sample.name, output_directory=str(sample.output), status='running', started=now(),
                   pipeline_log=str(logs / (sample.name + '.pipeline.log')),
                   publish_log=str(logs / (sample.name + '.publish.log')))
        rows.append(row)
        print('\n[' + str(index) + '/' + str(len(samples)) + '] ' + sample.name, flush=True)
        try:
            row['stage'] = 'verify' if args.publish_only else 'pipeline'
            save_summary(summary, rows)
            if not args.publish_only:
                # Recheck immediately before this sample, since earlier ones may take hours.
                audit_samples(args, [sample], enforce=True)
                run_logged(pipeline_command(args, sample), Path(row['pipeline_log']))
            row['stage'] = 'verify'
            save_summary(summary, rows)
            verify_completed(sample.output)
            row['stage'] = 'publish'
            save_summary(summary, rows)
            run_logged(publish_command(args, sample), Path(row['publish_log']))
            row.update(status='published', detail='s3://' + args.bucket + '/' + args.prefix.rstrip('/') + '/' + sample.output.name + '/')
            print('  Published: ' + sample.output.name, flush=True)
        except KeyboardInterrupt:
            row.update(status='interrupted', detail='Interrupted; rerun the same command to reuse verified stages.')
            raise
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
            row.update(status='failed', detail=str(exc))
            print('  FAILED: ' + str(exc), file=sys.stderr, flush=True)
        finally:
            row['finished'] = now()
            save_summary(summary, rows)
    failed = [row['sample'] for row in rows if row['status'] != 'published']
    print('\nPublished ' + str(len(rows) - len(failed)) + '/' + str(len(rows)) + ' samples. Summary: ' + str(summary), flush=True)
    return 1 if failed else 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--source', type=Path, default=Path.cwd(), help='Folder containing paired FASTQs and sample outputs (default cwd)')
    parser.add_argument('--pipeline-dir', type=Path, default=Path.home() / 'bin/MetaHopper')
    parser.add_argument('--publisher', type=Path, default=Path.home() / 'bin/metahopper-live/scripts/publish_reports.py')
    parser.add_argument('-c', '--config', type=Path, help='Default: PIPELINE_DIR/config.json')
    parser.add_argument('-t', '--threads', type=int, default=24)
    parser.add_argument('--memory', type=int, default=64, help='GB, one sample at a time')
    parser.add_argument('--pathway-jobs', type=int, default=3, help='Concurrent bins within the shared thread budget')
    parser.add_argument('--bucket', default=BUCKET)
    parser.add_argument('--prefix', default='reports/aphrophora/')
    parser.add_argument('--aws-profile')
    parser.add_argument('--graph-max-nodes', type=int, default=4000000)
    parser.add_argument('--megahit-max-links', type=int, default=10000000)
    parser.add_argument('--sample', action='append', help='Only this input stem; repeat to select several. Explicit APHCOR uses its existing -2 output.')
    parser.add_argument('--exclude', action='append', default=[], help='Additional input stem to exclude')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--plan', action='store_true', help='Print selected commands; no analysis, checks or uploads')
    mode.add_argument('--preflight-only', action='store_true', help='Check software, database, HMM panel and AWS credentials; no analysis/upload')
    mode.add_argument('--publish-only', action='store_true', help='Verify completed outputs and retry publishing without running the pipeline')
    mode.add_argument('--reuse-audit', action='store_true', help='Read saved settings/checkpoints and print reuse decisions; no analysis, AWS calls or output changes')
    parser.add_argument('--allow-recompute', action='store_true', help='Allow retained assemblies/searches to be rebuilt when their reuse proofs fail')
    parser.add_argument('--force-overwrite', action='store_true', help='Explicitly recompute enabled pipeline stages instead of reusing successes')
    args = parser.parse_args(argv)
    for key in ('source', 'pipeline_dir', 'publisher'):
        setattr(args, key, getattr(args, key).expanduser().resolve())
    args.config = (args.config or args.pipeline_dir / 'config.json').expanduser().resolve()
    if not args.source.is_dir():
        parser.error('--source must be an existing directory')
    if args.threads < 1 or args.memory < 4 or not 1 <= args.pathway_jobs <= args.threads:
        parser.error('Use threads >= jobs >= 1 and memory >= 4 GB')
    if args.graph_max_nodes < 1 or args.megahit_max_links < 1:
        parser.error('Graph limits must be positive')
    if args.publish_only and args.force_overwrite:
        parser.error('--publish-only cannot be combined with --force-overwrite')
    if not re.fullmatch(r'[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]', args.bucket):
        parser.error('--bucket must be an S3 bucket name')
    if (not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._/-]*/?', args.prefix)
            or any(p in ('', '.', '..') for p in args.prefix.rstrip('/').split('/'))):
        parser.error('Invalid S3 prefix')
    return args


def main(argv=None):
    args = parse_args(argv)
    excluded = set(args.exclude) | (set() if args.sample else {'APHCOR_AUC2351'})
    if args.sample and set(args.sample) & set(args.exclude):
        raise ValueError('A sample cannot be both selected and excluded.')
    samples = discover_samples(args.source, args.sample, excluded)
    print('Selected ' + str(len(samples)) + ' samples: ' + ', '.join(s.name for s in samples), flush=True)
    if args.plan:
        for sample in samples:
            print('\n' + shlex.join(pipeline_command(args, sample)))
            print(shlex.join(publish_command(args, sample)))
        return 0
    if args.reuse_audit:
        return 1 if audit_samples(args, samples) else 0
    with batch_lock(args.source):
        return run_batch(args, samples)


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print('\nBatch interrupted. Completed samples remain saved; rerun to resume.', file=sys.stderr)
        raise SystemExit(130)
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or str(exc)).strip()
        print('ERROR: ' + detail, file=sys.stderr)
        raise SystemExit(1)
    except (OSError, ValueError, RuntimeError, ImportError, subprocess.SubprocessError) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        raise SystemExit(1)
