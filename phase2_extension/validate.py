"""Standard-library validation; all artifacts confined to this extension."""
import sys
sys.dont_write_bytecode=True
import io
from pathlib import Path
import subprocess
import time
import unittest
import uuid
import pipeline as P
import scoring as S
import screen
import storage as A
H=P.H

def family():
    aa='ACDEFGHIKLMNPQRSTVWY'
    return {'rows':{n:(aa[k:]+aa[:k])*7 for k,n in enumerate('abcd')},'core_mask':[1]*140,
        'focal':['a','b'],'anchor_0based':35,'component':'toy','reference_set':'toy'}

def fake(cmd):
    rows=S.parse_fasta(Path(cmd[-1]).read_text());width=max(map(len,rows.values()))
    return P.B.fasta({n:s.ljust(width,'-') for n,s in rows.items()}).encode(),b'Converged.\n',0

class Validation(unittest.TestCase):
    def setUp(self):
        self.old=P.H;self.root=H/'test_artifacts'/uuid.uuid4().hex;self.root.mkdir(parents=True)
        for n in ['pipeline.py','scoring.py']:A.atomic(self.root/n,(H/n).read_bytes())
        P.H=self.root;P.B.HERE=self.root
    def tearDown(self):P.H=self.old;P.B.HERE=self.old
    def case(self,a='partial',fam=None):return P.B.prepare_case('toy',fam or family(),a)
    def test_deletion_lengths_overlap_union_and_coordinates(self):
        for a,o in [('none',0),('partial',6),('high',11)]:
            ref,inp,m=P.dataset(family(),a);self.assertEqual(m['shared_deleted_columns'],o);self.assertEqual(m['union_deleted_columns'],24-o)
            self.assertEqual(len(inp['a']),128);self.assertEqual(len(inp['b']),128);self.assertEqual(len(inp['c']),140)
    def test_survivor_mapping_exact_identity(self):
        f=family()
        for a in P.C['offsets']:
            ref,inp,m=P.dataset(f,a);S.validate_alignment(ref,inp)
            for n,coords in m['ancestry'].items():self.assertEqual(inp[n],''.join(f['rows'][n][c-1] for c in coords))
    def test_ambiguous_symbols_preserved_not_imputed(self):
        f=family();f['rows']['c']='XBZ'+f['rows']['c'][3:]
        root=self.case(fam=f);ref=S.parse_fasta((root/'reference.fasta').read_text());truth,b,known,m=P.read_case(root)
        self.assertTrue(ref['c'].startswith('XBZ'));self.assertTrue(all(not any(r in [('c',1),('c',2),('c',3)] for r in p) for p in truth))
        self.assertFalse(any(('c',i) in known for i in [1,2,3]));P.B.align(root,0,fake);P.score(root,0)
    def test_unknown_focal_deletion_is_rejected(self):
        f=family();f['rows']['a']=f['rows']['a'][:35]+'X'+f['rows']['a'][36:]
        with self.assertRaises(AssertionError):P.dataset(f,'none')
    def test_gap_focal_deletion_is_rejected(self):
        f=family();f['rows']['a']=f['rows']['a'][:35]+'-'+f['rows']['a'][36:]
        with self.assertRaises(AssertionError):P.dataset(f,'none')
    def test_screening_all_lengths_and_overlap_fractions(self):
        f=family()
        for d in [6,8,10,12]:
            selected=screen.site(f,d,10);self.assertIsNotNone(selected)
            for offset,overlap in [(d,0),(d//2,d//2),(1,d-1)]:
                self.assertEqual(len(set(range(d))&set(range(offset,offset+d))),overlap)
    def test_short_flanks_still_require_three_core_residues(self):
        f=family();s=screen.site(f,12,5);self.assertIsNotNone(s)
        f['core_mask']=[0]*140;self.assertIsNone(screen.site(f,12,5))
    def test_targeted_accepts_interrupted_noncore_flanks(self):
        f=family();f['core_mask']=[0]*140
        # 24 reliable deletion columns, with three reliable residues on each flank.
        for i in list(range(35,59))+list(range(30,33))+list(range(60,63)):f['core_mask'][i]=1
        self.assertIsNotNone(screen.site(f,12,10,'targeted'))
        self.assertIsNone(screen.site(f,12,10,'continuous'))
    def test_ambiguity_source_policy_increases_eligibility_without_letters_changed(self):
        f=family();f['rows']['d']='XBZ'+f['rows']['d'][3:]
        self.assertIsNone(screen.site(f,12,10,symbols='strict'))
        self.assertIsNotNone(screen.site(f,12,10,symbols='retained'))
    def test_shared_protein_components_include_excluded_bridges(self):
        d={'a':{'rows':{'x':'AC','y':'DE'}},'b':{'rows':{'y':'DE','z':'XB'}},'c':{'rows':{'w':'XB'}}}
        comp=screen.components(d);self.assertEqual(len(set(comp.values())),1)
    def test_hand_worked_scoring_and_inverse(self):
        truth=S.pairs({'a':'AC','b':'AC','c':'AC'});base=S.pairs({'a':'AC-','b':'-AC','c':'AC-'})
        r=S.compare(base,truth,truth,{('a',1),('a',2)});self.assertEqual(r['global']['recovered'],4);self.assertEqual(r['global']['lost'],0)
        self.assertEqual(r['local']['recovered'],2);self.assertEqual(S.compare(truth,base,truth,set())['global']['lost'],4)
    def test_repeat_letters_not_homology(self):self.assertEqual(len(S.pairs({'a':'AA-','b':'-AA'})&S.pairs({'a':'AA','b':'AA'})),0)
    def test_row_order_and_allgap_invariance(self):self.assertEqual(S.pairs({'a':'A-C','b':'A-C'}),S.pairs({'b':'AC','a':'AC'}))
    def test_local_and_global_bookkeeping(self):
        t=S.pairs({'a':'AC','b':'AC','c':'AC'});b=S.pairs({'a':'AC-','b':'-AC','c':'AC-'})
        r=S.compare(b,t,t,{('a',1)})
        for k in ['recovered','lost','reference_denominator']:self.assertEqual(r['global'][k],r['local'][k]+r['far'][k])
        self.assertIsNone(S.ratio(0,0))
    def test_only_core_standard_endpoints_scored(self):
        f=family();f['core_mask']=[0]*140
        for i in range(35,59):f['core_mask'][i]=1
        root=self.case(fam=f);truth,b,known,m=P.read_case(root)
        self.assertTrue(all(m['core_mask'][m['ancestry'][r[0]][r[1]-1]-1] for pair in truth for r in pair))
    def test_scoring_failure_preserves_raw_alignment(self):
        root=self.case();P.B.align(root,0,fake);p=root/'alignments/cap0/aligned.fasta';before=(A.sha(p),p.stat().st_mtime_ns)
        with self.assertRaises(RuntimeError):P.score(root,0,fail=True)
        self.assertFalse(P.B.align(root,0,fake));P.score(root,0);self.assertEqual(before,(A.sha(p),p.stat().st_mtime_ns))
    def test_corrupt_hash_and_marker_rerun_only_affected_alignment(self):
        root=self.case();P.B.align(root,0,fake);P.B.align(root,2,fake)
        A.atomic(root/'alignments/cap2/aligned.fasta','>a\nA\n');self.assertFalse(P.B.align(root,0,fake));self.assertTrue(P.B.align(root,2,fake))
        A.atomic(root/'alignments/cap2/done.json','{');self.assertTrue(P.B.align(root,2,fake))
    def test_partial_quarantine_and_invalid_output(self):
        root=self.case();A.atomic(root/'alignments/cap0/abandoned.partial','bad')
        with self.assertRaises(ValueError):P.B.align(root,0,lambda cmd:(b'>a\nA\n',b'bad',0))
        self.assertFalse((root/'alignments/cap0/done.json').exists());self.assertTrue(P.B.align(root,0,fake))
        self.assertTrue(list((root/'alignments/cap0/quarantine').glob('*/*.partial')))
    def test_repeated_scoring_idempotent(self):
        root=self.case();P.B.align(root,0,fake);P.score(root,0);p=root/'scores/cap0.json';before=(A.sha(p),p.stat().st_mtime_ns)
        self.assertFalse(P.B.align(root,0,fake));P.score(root,0);self.assertEqual(before,(A.sha(p),p.stat().st_mtime_ns))
    def test_baseline_hash_dependency_invalidates_only_score(self):
        root=self.case();P.B.align(root,0,fake);P.B.align(root,2,fake);P.score(root,2)
        old=A.read_json(root/'scores/cap2_done.json')['dependencies']['baseline']
        p=root/'alignments/cap0/aligned.fasta';A.atomic(p,p.read_text()+'\n');P.score(root,2)
        self.assertNotEqual(old,A.read_json(root/'scores/cap2_done.json')['dependencies']['baseline'])
    def test_different_config_not_reused(self):
        root=self.case();deps=A.read_json(root/'dataset_done.json')['dependencies'];deps['config']='changed'
        self.assertFalse(A.complete(root/'dataset_done.json',deps))
    def test_actual_terminated_process_skips_two_completed_alignments(self):
        roots=[self.case(a) for a in ['none','partial','high']]
        code="import sys,os;from pathlib import Path;sys.path.insert(0,sys.argv[1]);import pipeline as P;from validate import fake;P.H=Path(sys.argv[2]);P.B.HERE=P.H;P.B.align(P.H/'cases/toy/none',0,fake);P.B.align(P.H/'cases/toy/partial',0,fake);os._exit(73)"
        # pipeline imports the frozen config from its source directory, not cwd.
        p=subprocess.run([sys.executable,'-B','-c',code,str(H),str(self.root)],cwd=self.root,capture_output=True);self.assertEqual(p.returncode,73,p.stderr)
        before={str(r):(A.sha(r/'alignments/cap0/aligned.fasta'),(r/'alignments/cap0/aligned.fasta').stat().st_mtime_ns) for r in roots[:2]}
        calls=[]
        def invoke(cmd):calls.append(cmd);return fake(cmd)
        for root in roots:P.B.align(root,0,invoke)
        self.assertEqual(len(calls),1)
        for r,(digest,stamp) in before.items():self.assertEqual((digest,stamp),(A.sha(Path(r)/'alignments/cap0/aligned.fasta'),(Path(r)/'alignments/cap0/aligned.fasta').stat().st_mtime_ns))

def main():
    output=io.StringIO();t=time.time();result=unittest.TextTestRunner(stream=output,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Validation))
    ident=uuid.uuid4().hex;A.atomic(H/'validation_logs'/f'{ident}.log',output.getvalue())
    obj={'passed':result.wasSuccessful(),'tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'runtime_seconds':time.time()-t,'utc_epoch':time.time(),'log':f'validation_logs/{ident}.log','real_mafft_calls':0}
    A.json_write(H/'validation_logs'/f'{ident}.json',obj);A.json_write(H/'validation_latest.json',obj);print(output.getvalue());print(obj)
    return 0 if result.wasSuccessful() else 1
if __name__=='__main__':raise SystemExit(main())
