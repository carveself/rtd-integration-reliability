"""Audit packaged result counts, source hashes, reference flags and summaries."""
from pathlib import Path
import csv,json,hashlib,re
ROOT=Path(__file__).resolve().parent

def read(path):return list(csv.DictReader(path.open()))

def main():
    for folder,n in [('results',24),('results_stress',12)]:
        p=ROOT/folder
        assert len(read(p/'references.csv'))==n
        assert len(read(p/'benchmark.csv'))==3*n
        assert len(read(p/'timings.csv'))==9*n
        assert len(json.loads((p/'cases.json').read_text()))==n
        refs=read(p/'references.csv')
        failed=[r for r in refs if r['reference_accepted']=='False']
        if folder=='results':assert not failed
        else:
            audits=read(p/'precision_audit.csv');assert len(failed)==len(audits)==4
            assert {(r['id'],r['sample']) for r in failed}=={(r['id'],r['sample']) for r in audits}
            assert all(r['accepted_posthoc_audit']=='True' for r in audits)
        provenance=json.loads((p/'provenance.json').read_text())
        protocol='protocol.json' if folder=='results' else 'protocol_stress.json'
        assert hashlib.sha256((ROOT/protocol).read_bytes()).hexdigest()==provenance['protocol_sha256']
        for name,digest in provenance['files'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,(folder,name)
    assert len(read(ROOT/'results_legacy/benchmark.csv'))==33
    assert len(read(ROOT/'results_legacy/sample_errors.csv'))==924
    summary=json.loads((ROOT/'results/summary.json').read_text())
    assert summary['legacy_308_posthoc']['uniform']['undercovered_items']==5
    for method in ['uniform','spectral','pole']:
        assert summary['main_24'][method]['passed_on_accepted']==24
        assert summary['adversarial_12_posthoc'][method]['passed_on_accepted']==12
    for path in ROOT.rglob('*.md'):
        if '.git' in path.parts:continue
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in target or target.startswith('#'):continue
            assert (path.parent/target.split('#')[0]).exists(),(path,target)
    print('PASS: result counts, audited references, numerical claims, recorded source hashes and local document links')

if __name__=='__main__':main()
