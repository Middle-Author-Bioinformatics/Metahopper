import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import metahopper_annotate as annotation
from repair_bakta_db import diagnose


class BaktaDatabaseTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.root=Path(self.tmp.name)
        self.info={'major':5,'minor':0,'date':'2024-01-01'}
        (self.root/'psc.dmnd').write_text('fixture')
        (self.root/'bakta.db').write_text('fixture')
        self.save()
    def tearDown(self):self.tmp.cleanup()
    def save(self):(self.root/'version.json').write_text(json.dumps(self.info))
    def test_missing_type_is_detected_before_annotation_or_tool_execution(self):
        with patch.object(annotation.shutil,'which') as which:
            with self.assertRaisesRegex(ValueError,'missing type'):
                annotation.validate_config({'bakta_db':str(self.root)})
            which.assert_not_called()
    def test_legacy_type_can_only_be_inferred_after_schema_and_file_validation(self):
        metadata,kind=diagnose(self.root,5,['bakta.db'])
        self.assertEqual(kind,'full');self.assertNotIn('type',metadata)
        self.assertNotIn('type',json.loads((self.root/'version.json').read_text()))
        with self.assertRaisesRegex(ValueError,'does not match'):
            diagnose(self.root,6,['bakta.db'])
        with self.assertRaisesRegex(ValueError,'files missing'):
            diagnose(self.root,5,['missing.db'])
    def test_ambiguous_or_conflicting_type_is_rejected(self):
        (self.root/'pscc.dmnd').write_text('fixture')
        with self.assertRaisesRegex(ValueError,'unambiguously'):
            diagnose(self.root,5,['bakta.db'])
        (self.root/'pscc.dmnd').unlink()
        self.info['type']='light';self.save()
        with self.assertRaisesRegex(ValueError,'conflicts'):
            diagnose(self.root,5,['bakta.db'])
    def test_repair_cli_backs_up_metadata_and_preserves_permissions(self):
        package=self.root/'packages/bakta';package.mkdir(parents=True)
        (package/'__init__.py').write_text("__version__='1.10.3'\n__db_schema_version__=5\n")
        (package/'db.py').write_text("import json\nFILE_NAMES=['bakta.db']\ndef check(path):\n assert json.loads((path/'version.json').read_text())['type']=='full'\n")
        original=(self.root/'version.json').read_bytes();(self.root/'version.json').chmod(0o644)
        result=subprocess.run([sys.executable,str(Path(__file__).with_name('repair_bakta_db.py')),'--db',str(self.root),'--repair-missing-type'],
                              env={**os.environ,'PYTHONPATH':str(package.parent)},capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(json.loads((self.root/'version.json').read_text())['type'],'full')
        backups=list(self.root.glob('version.json.before-type-repair-*'))
        self.assertEqual(len(backups),1);self.assertEqual(backups[0].read_bytes(),original)
        self.assertEqual((self.root/'version.json').stat().st_mode&0o777,0o644)
    def test_incompatible_schema_never_changes_metadata(self):
        self.info['major']=4;self.save();original=(self.root/'version.json').read_bytes()
        with self.assertRaisesRegex(ValueError,'does not match'):
            diagnose(self.root,5,['bakta.db'])
        self.assertEqual((self.root/'version.json').read_bytes(),original)


if __name__=='__main__':unittest.main()
