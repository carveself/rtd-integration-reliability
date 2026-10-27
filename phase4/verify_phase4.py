from pathlib import Path
import csv,json,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parent

def read(p):return list(csv.DictReader(p.open()))

def main():
    assert len(read(P/'results/benchmark.csv'))==96
    assert len(read(P/'results/timings.csv'))==288
    refs=read(P/'results/references.csv');assert len(refs)==24
    assert all(r['accepted']=='True' for r in refs)
    assert sum(r['reference_type']=='hybrid_MP_repair' for r in refs)==3
    assert len(read(P/'results_regression/benchmark.csv'))==44
    assert len(read(P/'results_regression/sample_errors.csv'))==1232
    summary=json.loads((P/'results/summary.json').read_text())
    for m in ['uniform','spectral','pole','selective']:assert summary['fresh'][m]['passed']==24
    assert [summary['legacy'][m]['undercovered'] for m in ['uniform','spectral','pole','selective']]==[5,0,0,3]
    assert summary['fresh']['selective']['energy_evaluations']==60714
    assert summary['fresh']['selective']['energy_evaluations']<summary['fresh']['uniform']['energy_evaluations']
    for name,digest in json.loads((P/'results/provenance.json').read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    print('PASS: fresh/legacy result counts, all reference audits, 3 retained regression failures and recorded source hashes')

if __name__=='__main__':main()
