import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import metahopper_pathway_view as view
import metahopper_pathways as backend
import metahopper_report as report
from pathway_test_fixture import fixture


class PathwayTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)/'Sample'
        self.output=fixture(self.root)

    def tearDown(self):self.tmp.cleanup()

    def data(self):return view.read_results(self.output,self.root)

    def test_full_evidence_retains_zero_failed_and_stale_bins(self):
        data=self.data();self.assertEqual(len(data['summary']),4)
        summaries={r['bin_id']:r for r in data['summary']}
        self.assertEqual(summaries['A']['panel_percent'],100)
        self.assertEqual(summaries['A']['route_status'],'unresolved_unmeasured_steps')
        self.assertEqual(summaries['Zero']['panel_percent'],0)
        self.assertIsNone(summaries['Failed']['panel_percent'])
        self.assertEqual(summaries['Stale']['panel_status'],'stale_input')
        self.assertIsNone(summaries['Stale']['panel_percent'])
        self.assertFalse(any(h['bin_id']=='Stale' for h in data['hits']))

    def test_sequences_ids_shared_review_and_unmeasured_calls_are_preserved(self):
        data=self.data();self.assertEqual(len(data['hits']),3)
        self.assertEqual(sum(h['status']=='accepted' for h in data['hits']),2)
        self.assertEqual(sum(h['status']!='accepted' or h['target_call']!='specific' for h in data['hits']),2)
        self.assertEqual(data['proteins']['A\x1fc0000001_1']['aa'],'MK*')
        self.assertEqual(data['proteins']['A\x1fc0000001_1']['nt'],'ATGAAATAA')
        self.assertEqual(next(r for r in data['components'] if r['bin_id']=='A' and r['component']=='c3')['value'],'NA')
        self.assertEqual(next(r for r in data['bins'] if r['bin_id']=='A')['accepted_proteins'],1)

    def test_legacy_results_are_discovered_without_running_tools(self):
        legacy=self.root/'final/bins/pathwaygenie';shutil.move(str(self.output),legacy)
        layout=report.RunLayout(self.root,'final');data=view.collect_pathways(layout)
        self.assertEqual(data['directory'],'final/bins/pathwaygenie')
        fragment=report.render_annotation(layout)
        self.assertIn('id="annotation-pathways-tab"',fragment)
        self.assertIn('id="pathway-heatmap"',fragment)

    def test_unfinished_or_locked_results_are_not_presented_as_complete(self):
        (self.output/'.pathwaygenie.lock').write_text('123')
        with self.assertRaisesRegex(ValueError,'locked'):self.data()
        (self.output/'.pathwaygenie.lock').unlink()
        path=self.output/'run_manifest.json';m=json.loads(path.read_text());m.pop('finished_utc');path.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError,'finished'):self.data()

    def test_active_wrapper_is_locked_but_a_released_lock_is_safe(self):
        with backend.result_lock(self.output):
            with self.assertRaisesRegex(ValueError,'active'):self.data()
        self.assertEqual(len(self.data()['bins']),4)

    def test_different_rule_file_is_rejected(self):
        (self.output/'config_used.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'checksum'):self.data()

    def test_duplicate_summary_pairs_are_rejected(self):
        p=self.output/'pathway_summary.long.csv';lines=p.read_text().splitlines();p.write_text('\n'.join(lines+[lines[1]])+'\n')
        with self.assertRaisesRegex(ValueError,'duplicate'):self.data()

    def test_symlink_outside_sample_is_rejected(self):
        p=self.output/'all_model_hits.csv';outside=Path(self.tmp.name)/'outside.csv';p.replace(outside);p.symlink_to(outside)
        with self.assertRaisesRegex(ValueError,'outside'):self.data()

    def test_embedded_data_escapes_script_terminators_and_preserves_text(self):
        data=self.data();data['hits'][0]['source_product']='</script><img src=x onerror=alert(1)>'
        content=view.render_pathways(data);self.assertNotIn('</script><img',content)
        parsed=json.loads(re.search(r'id="metahopper-pathways-data">(.*?)</script>',content,re.S)[1])
        self.assertEqual(parsed['hits'][0]['source_product'],data['hits'][0]['source_product'])

    def test_asset_hashes_match_and_input_snapshot_is_sample_relative(self):
        data=self.data();self.assertTrue(data['files']);self.assertEqual(len(data['inputs']),4)
        for item in data['files']:
            p=self.root/item['path'];self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),item['sha256'])
            self.assertEqual(p.stat().st_size,item['size'])
        self.assertTrue(all(not Path(r['path']).is_absolute() for r in data['inputs']))

    def test_import_preserves_originals_and_requires_matching_sample(self):
        source=Path(self.tmp.name)/'outside-results';shutil.move(str(self.output),source)
        data=backend.import_results(self.root,source)
        self.assertEqual(data['directory'],'final/pathways/genus')
        self.assertTrue((source/'accepted_hits.csv').exists() or (source/'all_model_hits.csv').exists())
        with self.assertRaisesRegex(ValueError,'already exists'):backend.import_results(self.root,source)
        backend.import_results(self.root,source,force=True)
        self.assertEqual(len(list(self.output.parent.glob('genus.before-import-*'))),1)

    def test_import_wrong_sample_never_replaces_existing_output(self):
        source=Path(self.tmp.name)/'wrong-output';shutil.copytree(self.output,source)
        p=source/'run_manifest.json';m=json.loads(p.read_text());
        for item in m['bin_inputs']:item['sha256']='f'*64
        p.write_text(json.dumps(m));before=(self.output/'run_manifest.json').read_bytes()
        with self.assertRaisesRegex(ValueError,'No PathwayGenie input'):backend.import_results(self.root,source,force=True)
        self.assertEqual((self.output/'run_manifest.json').read_bytes(),before)

    def test_invalid_extra_arguments_cannot_override_managed_inputs(self):
        hmm=Path(self.tmp.name)/'HMMs';hmm.mkdir();(hmm/'download_manifest.json').write_text('{}')
        config=dict(pathwaygenie_hmm_dir=str(hmm),tools={'pathwaygenie':[sys.executable]},pathwaygenie_args=['--out=somewhere'])
        with self.assertRaisesRegex(ValueError,'override'):backend.runtime_config(config)

    def test_run_forwards_resources_force_and_retains_partial_evidence(self):
        hmm=Path(self.tmp.name)/'HMMs';hmm.mkdir();(hmm/'download_manifest.json').write_text('{}')
        observed=[]
        class Result:returncode=1
        def fake_run(command,**kwargs):observed.append(command);return Result()
        with patch.object(backend,'preflight_runtime',return_value=([sys.executable,'fake.py'],hmm,self.output/'config_used.json',[])),patch.object(backend.subprocess,'run',side_effect=fake_run):
            result=backend.run_pathways(self.root,threads=8,jobs=2,force=False,log=lambda _:None)
            self.assertEqual(result['status'],'partial');self.assertEqual(result['failed_bins'],['Failed'])
            self.assertNotIn('--force',observed[-1]);self.assertIn('8',observed[-1]);self.assertIn('2',observed[-1])
            backend.run_pathways(self.root,threads=8,jobs=2,force=True,log=lambda _:None)
            self.assertIn('--force',observed[-1])
        state=json.loads((self.root/'pathwaygenie.genus.status.json').read_text());self.assertEqual(state['status'],'partial')

    def test_preflight_checks_the_selected_python_environment_before_heavy_steps(self):
        hmm=Path(self.tmp.name)/'HMMs';hmm.mkdir();(hmm/'download_manifest.json').write_text('{}')
        config=dict(pathwaygenie_hmm_dir=str(hmm),pathwaygenie_env='')
        with patch.object(backend,'validate_panel',return_value={}),patch.object(backend.subprocess,'run',return_value=__import__('subprocess').CompletedProcess([],0,stdout='OK')) as run:
            backend.preflight_runtime(config,threads=4,jobs=2)
        self.assertEqual(run.call_count,2)
        self.assertIn('--version',run.call_args_list[0].args[0])
        self.assertIn('-c',run.call_args_list[1].args[0])
        self.assertIn('hmmsearch',run.call_args_list[1].args[0][-1])

    def test_missing_runtime_tools_fail_with_an_environment_message(self):
        hmm=Path(self.tmp.name)/'HMMs';hmm.mkdir();(hmm/'download_manifest.json').write_text('{}')
        config=dict(pathwaygenie_hmm_dir=str(hmm),pathwaygenie_env='')
        error=__import__('subprocess').CalledProcessError(1,['python'],output='Missing PathwayGenie tools: hmmsearch')
        with patch.object(backend,'validate_panel',return_value={}),patch.object(backend.subprocess,'run',side_effect=error):
            with self.assertRaisesRegex(ValueError,'environment/dependency.*hmmsearch'):
                backend.preflight_runtime(config)

    def test_default_install_uses_nested_resources_and_active_python(self):
        home=Path(self.tmp.name)/'relocated-repo';nested=home/'pathwaygenie';nested.mkdir(parents=True)
        (nested/'PathwayGenie.py').write_text('# fixture')
        (nested/'pathwaygenie_config.json').write_text('{}')
        hmm=nested/'Pathway_Gene_HMMs';hmm.mkdir();(hmm/'download_manifest.json').write_text('{}')
        with patch.object(backend,'HERE',home):
            prefix,found,rules,_=backend.runtime_config({},threads=4,jobs=2)
        self.assertEqual(prefix,[sys.executable,str(nested/'PathwayGenie.py')])
        self.assertEqual(found,hmm)
        self.assertEqual(rules,nested/'pathwaygenie_config.json')
        self.assertNotIn('conda',prefix)

    def test_bundled_panel_has_114_verified_profiles_and_three_libraries(self):
        prefix,hmm,rules,extra=backend.runtime_config({})
        result=backend.validate_panel(Path(prefix[-1]),hmm,rules,extra)
        self.assertEqual(result['profiles'],114)
        self.assertEqual(result['panel_profiles'],114)
        self.assertEqual(result['libraries'],3)

    def test_corrupted_bundled_source_hmm_fails_integrity_check(self):
        prefix,hmm,rules,_=backend.runtime_config({})
        damaged=Path(self.tmp.name)/'damaged-panel';shutil.copytree(hmm,damaged)
        manifest=json.loads((damaged/'download_manifest.json').read_text())
        source=damaged/manifest['models'][0]['relative_file']
        with source.open('ab') as handle:handle.write(b'\n# changed\n')
        with self.assertRaisesRegex(ValueError,'checksum'):
            backend.validate_panel(Path(prefix[-1]),damaged,rules)

    def test_permanent_report_contains_three_main_tabs_and_scoped_download_fix(self):
        args=__import__('argparse').Namespace(title=None,min_length=1000,max_points=100,include_uncoordinated=False,no_fasta_stats=False,skip_graphs=True)
        content=report.build_report(report.RunLayout(self.root,'final'),args)
        self.assertIn('"id": "annotation"',content)
        self.assertNotIn('"id": "pathways", "label"',content)
        self.assertNotIn('__DOWNLOAD_LAYOUT_CSS__',content)
        self.assertIn('flex-wrap:nowrap!important',content)
        self.assertIn('"source_state":"stale_input"',content)


if __name__=='__main__':unittest.main()
