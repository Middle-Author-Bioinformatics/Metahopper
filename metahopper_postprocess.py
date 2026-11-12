#!/usr/bin/env python3
"""Assess and annotate saved final bins without rebuilding them or remapping reads."""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import MetaHopper as pipeline
import metahopper_annotate as annotation
import metahopper_pathways as pathways
from metahopper_state import Checkpoints, fingerprint, output_lock

HERE = Path(__file__).resolve().parent


def load_config(path):
    path = Path(path).expanduser().resolve()
    config = json.loads(path.read_text())
    for key in ('diamond_db', 'taxdump', 'phyloflash_db', 'bakta_db', 'annotation_metadata',
                'pathwaygenie_script', 'pathwaygenie_hmm_dir', 'pathwaygenie_config'):
        if config.get(key):
            value = Path(os.path.expandvars(config[key])).expanduser()
            config[key] = str((path.parent / value).resolve())
    return config


def validate_final(root):
    root = Path(root).resolve()
    manifest = json.loads((root / 'run.json').read_text())
    if manifest.get('status') != 'completed':
        raise ValueError('Existing-final mode requires a completed saved run: ' + str(root))
    base = root / 'final'
    fastas = sorted((base / 'bins/genus').glob('*.fasta'))
    if not fastas or any(not p.stat().st_size for p in fastas):
        raise ValueError('Existing nonempty final genus FASTAs are required: ' + str(base / 'bins/genus'))
    for path in (base / 'bins').glob('*/*.fasta'):
        if not path.resolve().is_relative_to(root):
            raise ValueError('Final bin FASTA points outside the sample: ' + str(path))
    summaries = pipeline.read_summaries(base, ['genus'])
    if set(summaries['genus']) != {p.stem for p in fastas}:
        raise ValueError('Final genus summary and FASTA names do not agree: ' + str(base))
    for relative in ('classification/contig_classification.tsv', 'binarena/binarena_input.tsv'):
        if not (base / relative).is_file():
            raise ValueError('Saved final report input is missing: ' + str(base / relative))
    return fastas


def validate(config, threads=24, jobs=3):
    annotation.validate_config(config)
    pathways.preflight_runtime(config, threads, jobs)
    prefix = config.get('tools', {}).get('checkm', ['checkm'])
    if not isinstance(prefix, list) or not prefix or not all(isinstance(v, str) and v for v in prefix):
        raise ValueError('tools.checkm must be an executable-prefix array')
    result = subprocess.run(prefix + ['lineage_wf', '-h'], text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=60, check=True)
    if not all(flag in result.stdout for flag in ('--reduced_tree', '--pplacer_threads', '--tab_table')):
        raise ValueError('The configured CheckM must support the existing lineage_wf interface (as CheckM 1.2.3 does).')


def ensure_final_export(root, fastas):
    """Provide the report's canonical FASTA without changing any bin sequences."""
    canonical = root / 'final/assembly/consolidated_contigs.fasta'
    if canonical.is_file():
        return canonical
    canonical.parent.mkdir(parents=True, exist_ok=True)
    legacy = canonical.with_name('contigs.fasta')
    with tempfile.NamedTemporaryFile(dir=canonical.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        if legacy.is_file():
            shutil.copy2(legacy, temporary)
        else:
            seen = set()
            with temporary.open('w') as handle:
                for path in fastas:
                    for identifier, sequence in pipeline.fasta(path):
                        if identifier in seen:
                            raise ValueError('Duplicate contig ID across final genus bins: ' + identifier)
                        seen.add(identifier)
                        handle.write('>' + identifier + '\n' + sequence + '\n')
            if not seen:
                raise ValueError('Final genus FASTAs contain no sequences')
        temporary.replace(canonical)
    finally:
        temporary.unlink(missing_ok=True)
    return canonical


def sequence_snapshot(root):
    return {str(p.relative_to(root)): fingerprint(p) for p in sorted((root / 'final/bins').glob('*/*.fasta'))}


def prepare_report(root, max_nodes, max_links):
    graph_command = [sys.executable, str(HERE / 'metahopper_graph.py'), '-i', str(root),
                     '--bin-set', 'final', '--max-nodes', str(max_nodes)]
    if (root / 'megahit/final.contigs.fa').is_file() and not (root / 'megahit/assembly.gfa').is_file():
        graph_command += ['--export-megahit-graph', '--megahit-max-links', str(max_links)]
    subprocess.run(graph_command, check=True)
    subprocess.run([sys.executable, str(HERE / 'metahopper_report.py'), '-i', str(root),
                    '--bin-set', 'final', '--graph-max-nodes', str(max_nodes), '--pathways-rank', 'genus'], check=True)


def run(root, config_path, threads=24, jobs=3, memory=64, max_nodes=4000000, max_links=10000000):
    root = Path(root).resolve()
    config_path = Path(config_path).expanduser().resolve()
    config = load_config(config_path)
    with output_lock(root):
        fastas = validate_final(root)
        frozen = sequence_snapshot(root)
        validate(config, threads, jobs)
        ensure_final_export(root, fastas)
        options = pipeline.parse_args(['--resume', str(root), '-c', str(config_path),
                                      '-t', str(threads), '--memory', str(memory), '--checkm'])
        contexts = pipeline.checkpoint_context(options, config)
        summaries = pipeline.read_summaries(root / 'final', ['genus'])
        old_contexts = None
        previous = json.loads((root / 'run.json').read_text())
        try:
            old_contexts = pipeline.checkpoint_context(argparse.Namespace(**previous['arguments']), previous['config'])
        except (KeyError, TypeError, AttributeError):
            pass
        # A valid full-pipeline quality result can satisfy the genus assessment.
        main_checks = Checkpoints.read_only(root, contexts=contexts, previous_contexts=old_contexts)
        all_ranks = [rank for rank in options.ranks if (root / 'final/bins' / rank / 'summary.tsv').is_file()]
        table = root / 'final/checkm.tsv'
        inputs = [root / 'final/assembly/consolidated_contigs.fasta'] + [
            p for p in pipeline.bin_outputs(root / 'final', all_ranks) if p.suffix == '.fasta']
        existing_quality = main_checks.peek('checkm', inputs=inputs,
                                            required=[table] + [root / 'final/bins' / r / 'summary.tsv' for r in all_ranks],
                                            legacy_path=root / 'run.json')
        options.ranks = ['genus']
        contexts['checkm'] = dict(contexts['checkm'], ranks=['genus'])
        ck = Checkpoints(root, contexts=contexts, log=print, namespace='postprocess')
        ck.status.update(mode='existing_final_bins', bin_fastas=frozen)
        ck._save()
        warnings = []
        runner = pipeline.Runner(config, root, memory)
        try:
            print('Using saved final bin FASTAs; no QC, assembly, preliminary NR DIAMOND, mapping or Unicycler jobs.', flush=True)
            if existing_quality:
                ck.start('checkm', 'Reusing verified existing CheckM results')
                ck.done('checkm', [table, root / 'final/bins/genus/summary.tsv'])
            else:
                pipeline.run_checkm(runner, options, root, root / 'final/bins', ['genus'], summaries, warnings, ck)
            if ck.status['steps'].get('checkm', {}).get('status') == 'failed':
                raise RuntimeError('CheckM failed; inspect the sample CheckM log')
            ck.start('annotation', 'Bakta checks per-bin content caches')
            result = annotation.annotate_run(root, config, rank='genus', stage='final', threads=threads,
                                             mode=options.annotation_mode, timeout_h=options.annotation_timeout_h,
                                             translation_table=options.annotation_translation_table, log=print)
            if any(row['status'] == 'failed' for row in result['bins']):
                raise RuntimeError('Bakta has failed bins; successful bin annotations remain cached')
            ck.done('annotation', [result['summary']])
            ck.start('pathways', 'PathwayGenie checks per-bin content caches')
            result = pathways.run_pathways(root, config, rank='genus', stage='final', threads=threads, jobs=jobs, log=print)
            if result['status'] != 'completed':
                raise RuntimeError('PathwayGenie has failed bins; successful searches remain cached')
            ck.done('pathways', [result['manifest']])
            if sequence_snapshot(root) != frozen:
                raise RuntimeError('Final bin FASTAs changed during processing; publication is blocked')
            ck.start('report')
            prepare_report(root, max_nodes, max_links)
            ck.done('report', [root / 'metahopper_report.html'])
            if sequence_snapshot(root) != frozen:
                raise RuntimeError('Final bin FASTAs changed during report generation')
            ck.finish('completed')
        except BaseException as exc:
            ck.fail(error=str(exc) or type(exc).__name__)
            ck.finish('failed')
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', type=Path, required=True)
    parser.add_argument('-c', '--config', type=Path, required=True)
    parser.add_argument('-t', '--threads', type=int, default=24)
    parser.add_argument('--memory', type=int, default=64)
    parser.add_argument('--pathway-jobs', type=int, default=3)
    parser.add_argument('--graph-max-nodes', type=int, default=4000000)
    parser.add_argument('--megahit-max-links', type=int, default=10000000)
    args = parser.parse_args()
    if args.memory < 4 or not 1 <= args.pathway_jobs <= args.threads or min(args.graph_max_nodes, args.megahit_max_links) < 1:
        parser.error('Invalid thread, memory or graph limits')
    try:
        run(args.input, args.config, args.threads, args.pathway_jobs, args.memory, args.graph_max_nodes, args.megahit_max_links)
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        parser.exit(1, 'ERROR: ' + str(exc) + '\n')


if __name__ == '__main__':
    main()
