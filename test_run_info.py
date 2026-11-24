import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from metahopper_run_info import run_files


class RunInfoTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / 'Sample'
        self.root.mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, name, data):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def test_settings_and_all_saved_log_bytes_are_preserved_without_mutating_inputs(self):
        settings = dict(arguments={'checkm': True, 'threads': 24}, config={'tools': {'checkm': ['conda', 'run', '-n', 'checkm_env', 'checkm']}})
        run = self.write('run.json', json.dumps(settings).encode())
        full = (b'Complete stdout and stderr\n' * 50000) + b'END OF TOOL LOG\n'
        paths = [run, self.write('pipeline.log', b'Pipeline timeline\n'),
                 self.write('logs/tool.log', full), self.write('logs/commands.jsonl', b'{"command":["tool","--flag"]}\n'),
                 self.write('final/pathways/genus/bins/Test/hmm_main.command.log', b'HMM execution\n')]
        self.write('final/pathways/genus/config_used.json', b'{"minimum_coverage":0.6}')
        original = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in paths}
        links = run_files(self.root, 'final')
        info = json.loads((self.root / links['info']).read_text())
        output = (self.root / links['log']).read_bytes()
        self.assertEqual(info['settings']['run.json']['data'], settings)
        self.assertEqual(info['source_run'], 'Sample')
        self.assertEqual(len(info['logs']), 4)
        self.assertIn(full, output)
        self.assertIn(b'HMM execution', output)
        self.assertIn(b'"--flag"', output)
        for row in info['logs']:
            source = (self.root / row['path']).read_bytes()
            self.assertEqual(row['sha256'], hashlib.sha256(source).hexdigest())
            self.assertEqual(row['size'], len(source))
        for path, snapshot in original.items():
            self.assertEqual((path.read_bytes(), path.stat().st_mtime_ns), snapshot)
        before = (self.root / links['info']).read_bytes()
        run_files(self.root, 'final')
        self.assertEqual((self.root / links['info']).read_bytes(), before)

    def test_missing_records_and_logs_are_not_invented(self):
        links = run_files(self.root, 'preliminary')
        self.assertIsNone(links['log'])
        info = json.loads((self.root / links['info']).read_text())
        self.assertEqual(info['settings'], {})
        self.assertEqual(info['logs'], [])
        self.assertIn('No saved execution logs', (self.root / 'report_support/full_run.log').read_text())

    def test_external_metadata_logs_and_support_directory_are_rejected(self):
        outside = Path(self.tmp.name) / 'outside'
        outside.mkdir()
        secret = outside / 'secret.log'
        secret.write_text('must not be published')
        for name in ('run.json', 'logs/leak.log'):
            path = self.root / name
            path.parent.mkdir(exist_ok=True)
            path.symlink_to(secret)
            with self.assertRaisesRegex(ValueError, 'outside this run'):
                run_files(self.root, 'final')
            path.unlink()
        import shutil
        shutil.rmtree(self.root / 'report_support')
        (self.root / 'report_support').symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'outside this run'):
            run_files(self.root, 'final')


if __name__ == '__main__':
    unittest.main()
