#!/usr/bin/env python3
"""Diagnose a Bakta database; optionally repair only a demonstrably missing type."""
import argparse
import json
import os
import stat
from pathlib import Path
import tempfile
import time


def diagnose(db, expected_schema, required_files):
    info = json.loads((db / 'version.json').read_text())
    if not isinstance(info, dict):
        raise ValueError('version.json must contain a JSON object')
    if type(info.get('major')) is not int or info['major'] != expected_schema:
        raise ValueError(f'Database schema {info.get("major")!r} does not match installed Bakta schema {expected_schema}. '
                         'Download a compatible database into a new directory with bakta_db download --type full --output NEW_DIRECTORY')
    if not isinstance(info.get('date'), str) or type(info.get('minor')) is not int:
        raise ValueError('Missing/invalid date or minor version; restore authoritative metadata or download a compatible database')
    present = {kind for kind, filename in [('full', 'psc.dmnd'), ('light', 'pscc.dmnd')] if (db / filename).is_file()}
    if len(present) != 1:
        raise ValueError('Database type cannot be established unambiguously from psc.dmnd / pscc.dmnd')
    kind = next(iter(present))
    if info.get('type', kind) != kind:
        raise ValueError('Declared database type conflicts with its sequence database files')
    missing = [name for name in required_files if not (db / name).is_file()]
    if missing:
        raise ValueError('Database files missing: ' + ', '.join(missing))
    return info, kind


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', type=Path, required=True)
    parser.add_argument('--repair-missing-type', action='store_true', help='Back up and add type only after schema and file checks pass')
    args = parser.parse_args()
    try:
        import bakta
        import bakta.db as database
        db = args.db.resolve()
        info, kind = diagnose(db, bakta.__db_schema_version__, database.FILE_NAMES)
        print(json.dumps({'bakta': bakta.__version__, 'required_schema': bakta.__db_schema_version__,
                          'database': str(db), 'metadata': info, 'detected_type': kind}, indent=2))
        if 'type' in info:
            database.check(db)
            print('Database metadata and required files validated. No changes made.')
            return
        if not args.repair_missing_type:
            print('Schema and required files match. Re-run with --repair-missing-type to add the missing type with a backup.')
            return
        original = (db / 'version.json').read_bytes()
        backup = db / f'version.json.before-type-repair-{time.time_ns()}'
        backup.write_bytes(original)
        info['type'] = kind
        with tempfile.NamedTemporaryFile(mode='w', dir=db, prefix='version.', suffix='.tmp', delete=False) as handle:
            json.dump(info, handle, indent=2); handle.write('\n'); temporary = Path(handle.name)
        temporary.chmod(stat.S_IMODE((db / 'version.json').stat().st_mode))
        try:
            os.replace(temporary, db / 'version.json')
            try:
                database.check(db)
            except BaseException:
                (db / 'version.json').write_bytes(original)
                raise
        finally:
            temporary.unlink(missing_ok=True)
        print(f'Added verified type={kind}. Original metadata backed up at {backup}. Database sequence files were unchanged.')
    except (ImportError, OSError, ValueError) as exc:
        parser.exit(1, f'ERROR: {exc}\n')


if __name__ == '__main__':
    main()
