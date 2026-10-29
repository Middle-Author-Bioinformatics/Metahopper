"""Small, atomic completion records and a human-readable run log (standard library)."""
from __future__ import annotations
import hashlib
import json
import os
import time
from contextlib import contextmanager
from pathlib import Path


def fingerprint(path, recursive=False):
    """Fast file identity: path, size and nanosecond mtime; no scan of large reads/DBs."""
    if not path:
        return None
    path = Path(path).resolve()
    result = {'path': str(path)}
    try:
        stat = path.stat()
    except OSError:
        return dict(result, missing=True)
    result.update(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    if path.is_dir():
        candidates = sorted(p for p in path.rglob('*') if p.is_file()) if recursive else [
            path / name for name in ('version.json', 'db.json', 'nodes.dmp', 'names.dmp', 'merged.dmp')
            if (path / name).is_file()]
        result['files'] = [fingerprint(p) for p in candidates]
    return result


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def sequence_digest(records):
    value = hashlib.sha256()
    for name, sequence in records:
        value.update(str(name).encode()); value.update(b'\0')
        value.update(sequence.encode()); value.update(b'\0')
    return value.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, default=str) + '\n')
    temporary.replace(path)


def write_if_changed(path, text):
    """Keep stable metadata when regenerated bookkeeping has identical contents."""
    path = Path(path)
    if path.is_file() and path.read_text() == text:
        return
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(text)
    temporary.replace(path)


@contextmanager
def output_lock(out):
    """Prevent two pipeline processes from writing the same output directory."""
    import fcntl
    with (Path(out) / 'pipeline.lock').open('a+') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'Another MetaHopper process is using {out}') from None
        try:
            handle.seek(0); handle.truncate()
            handle.write(f'pid={os.getpid()}\n'); handle.flush()
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


class Checkpoints:
    """Only reuse explicit successful completion plus matching settings and outputs.

    The root JSON records current/completed steps; pipeline.log is an append-only
    text timeline. Per-step checkpoints retain output metadata for quick checks.
    Interrupted, failed, missing, changed, or corrupt records are never successes.
    """

    def __init__(self, out, enabled=True, force=False, contexts=None, previous_contexts=None, log=print):
        self.out = Path(out).resolve()
        self.dir = self.out / 'work/checkpoints'
        self.dir.mkdir(parents=True, exist_ok=True)
        self.enabled, self.force = enabled, force
        self.contexts, self.previous_contexts = contexts or {}, previous_contexts
        self.log = log
        self.status_path = self.out / 'pipeline_status.json'
        self.timeline = self.out / 'pipeline.log'
        self.legacy_run = not self.status_path.exists()
        try:
            self.status = json.loads(self.status_path.read_text())
            if not isinstance(self.status.get('steps'), dict):
                raise ValueError('Invalid status')
        except (OSError, ValueError, TypeError, AttributeError):
            self.status = {'schema': 1, 'steps': {}}
        previous = self.status.get('current_step')
        if previous:
            self._event(previous, 'interrupted', 'Previous process stopped before recording completion')
        self.status.update(status='running', current_step=None, started_at=self.now())
        self.active_keys = {}
        self._save()

    @staticmethod
    def now():
        return time.strftime('%FT%T%z')

    def _save(self):
        self.status['updated_at'] = self.now()
        atomic_json(self.status_path, self.status)

    def _event(self, stage, state, message='', **details):
        record = self.status['steps'].get(stage, {})
        record.update(status=state, updated_at=self.now(), message=message, **details)
        self.status['steps'][stage] = record
        with self.timeline.open('a') as handle:
            handle.write(f'{self.now()}\t{stage}\t{state}\t{message.replace(chr(10), " ")}\n')

    def _context(self, stage, previous=False):
        contexts = self.previous_contexts if previous else self.contexts
        if contexts is None:
            return None
        return contexts.get(stage, contexts.get(stage.split('__', 1)[0], {}))

    def _key(self, stage, inputs, options):
        return digest({'cache_schema': 1, 'context': self._context(stage),
                       'inputs': [fingerprint(p) for p in inputs if p], 'options': options})

    def _path(self, stage):
        if '/' in stage or '\\' in stage or stage in ('.', '..'):
            raise ValueError('Invalid checkpoint name')
        return self.dir / f'{stage}.json'

    def _lookup(self, stage, inputs=(), options=None, required=(), legacy_path=None):
        key = self._key(stage, inputs, options)
        if not self.enabled or self.force:
            return None, key
        marker = self._path(stage)
        if not marker.is_file() and legacy_path and self.legacy_run:
            marker = Path(legacy_path)
        try:
            record = json.loads(marker.read_text())
            if not isinstance(record, dict):
                return None, key
            if record.get('status') not in (None, 'completed'):
                return None, key
            if record.get('outcome', '').startswith('failed'):
                return None, key
            if stage == 'checkm' and any('CheckM failed' in str(message) for message in record.get('warnings', [])):
                return None, key
            if stage == 'inspection' and (record.get('errors') or any(
                    row.get('status') == 'failed' for row in record.get('trees', []))):
                return None, key
            paths = list(record.get('files', [])) + list(required)
            if not all(Path(p).is_file() for p in paths):
                return None, key
            if record.get('signature'):
                if record['signature'] != key:
                    return None, key
                if any(fingerprint(item['path']) != item for item in record.get('output_stats', [])):
                    return None, key
            else:
                # Older versions wrote success markers after zero-exit commands.
                # Adopt only with compatible saved arguments and older inputs.
                if self.previous_contexts is None or self._context(stage, True) != self._context(stage):
                    return None, key
                cutoff = marker.stat().st_mtime_ns
                def newer(value):
                    if isinstance(value, dict):
                        return value.get('mtime_ns', 0) > cutoff or value.get('missing', False) or any(newer(v) for v in value.values())
                    if isinstance(value, (list, tuple)):
                        return any(newer(v) for v in value)
                    return False
                if newer(self._context(stage)) or any(newer(fingerprint(p)) for p in inputs if p):
                    return None, key
                record = dict(record, legacy_adopted=True)
            return record, key
        except (OSError, ValueError, TypeError, KeyError, AttributeError):
            return None, key

    def peek(self, stage, **kwargs):
        return self._lookup(stage, **kwargs)[0]

    def get(self, stage, inputs=(), options=None, required=(), legacy_path=None, mutable_files=()):
        record, key = self._lookup(stage, inputs, options, required, legacy_path)
        self.active_keys[stage] = key
        if record is not None:
            if record.get('legacy_adopted'):
                # Convert old proofs once; subsequent checks include output metadata.
                record = dict(record, schema=2, status='completed', signature=key,
                              files=list(dict.fromkeys(map(str, list(record.get('files', [])) + list(required)))))
                record['output_stats'] = [fingerprint(p) for p in record['files'] if str(p) not in set(map(str, mutable_files))]
                atomic_json(self._path(stage), record)
            self._event(stage, 'reused', 'Successful completion and outputs verified')
            self._save()
            self.log(f'     reuse: {stage}')
            return record
        self.start(stage, key=key)
        return None

    def start(self, stage, message='', key=None):
        if key is not None:
            self.active_keys[stage] = key
        self._path(stage).unlink(missing_ok=True)
        self.status['current_step'] = stage
        self._event(stage, 'running', message, started_at=self.now())
        self._save()

    def done(self, stage, files=(), mutable_files=(), **info):
        files, mutable_files = list(map(str, files)), list(map(str, mutable_files))
        if not all(Path(p).is_file() for p in files + mutable_files):
            raise RuntimeError(f'{stage} did not produce all required outputs')
        record = dict(info, schema=2, status='completed', signature=self.active_keys.get(stage),
                      files=files + mutable_files, output_stats=[fingerprint(p) for p in files], time=self.now())
        atomic_json(self._path(stage), record)
        self._event(stage, 'completed', completed_at=self.now())
        self.status['current_step'] = None
        self._save()

    def fail(self, stage=None, error='Failed'):
        stage = stage or self.status.get('current_step')
        if stage:
            self._path(stage).unlink(missing_ok=True)
            self._event(stage, 'failed', str(error))
        self.status['current_step'] = None
        self._save()

    def skip(self, stage, reason):
        self._event(stage, 'skipped', reason)
        self.status['current_step'] = None
        self._save()

    def finish(self, status='completed'):
        self.status.update(status=status, current_step=None, finished_at=self.now())
        self._save()
