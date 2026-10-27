"""Recompute summaries from raw data and render standalone scientific figures."""
from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'vendor')]
from resonance import make_device,discover
from negf import transport

METHODS=['uniform','spectral','pole']
LABELS={'uniform':'Uniform paired Gauss','spectral':'Spectral fixed scale','pole':'Spectral + pole width'}
COLORS={'uniform':'#757575','spectral':'#0072B2','pole':'#D55E00'}

def read(p):return list(csv.DictReader(p.open()))

def data(folder,repair=False):
    refs={(r['id'],r['sample']):r for r in read(folder/'references.csv')}
    audits={(r['id'],r['sample']):r for r in read(folder/'precision_audit.csv')} if repair else {}
    rows=read(folder/'benchmark.csv')
    for r in rows:
        key=r['id'],r['sample'];ref=refs[key]
        truth=float(ref['reference_A_m2']);accepted=ref['reference_accepted']=='True'
        if key in audits:
            truth=float(audits[key]['reference_MP_A_m2']);accepted=audits[key]['accepted_posthoc_audit']=='True'
        r.update(truth=truth,accepted=accepted,error=abs(float(r['J_A_m2'])-truth),
                 energy=int(r['energy_evaluations']),seconds=float(r['seconds_median']))
        r['rel']=r['error']/abs(truth);r['covered']=r['error']<=float(r['estimated_error_A_m2'])
        r['pass']=r['rel']<=.001+1e-6/abs(truth)
    summary={}
    for m in METHODS:
        rs=[r for r in rows if r['method']==m]
        summary[m]={'cases':len(rs),'accepted_references':sum(r['accepted'] for r in rs),
                    'passed_on_accepted':sum(r['accepted'] and r['pass'] for r in rs),
                    'covered_on_accepted':sum(r['accepted'] and r['covered'] for r in rs),
                    'max_relative_error':max(r['rel'] for r in rs),
                    'maximum_error_to_indicator':max(r['error']/float(r['estimated_error_A_m2']) for r in rs),
                    'energy_evaluations':sum(r['energy'] for r in rs),
                    'sum_per_case_median_seconds':sum(r['seconds'] for r in rs),
                    'pole_residual_calls':sum(int(r['pole_residual_calls']) for r in rs),
                    'sum_first_run_setup_seconds':sum(float(r['setup_seconds']) for r in rs)}
    return rows,summary


def main():
    main_rows,sm=data(ROOT/'results');stress_rows,ss=data(ROOT/'results_stress',True)
    legacy=read(ROOT/'results_legacy/benchmark.csv');sl={}
    for m in METHODS:
        a=[r for r in legacy if r['method']==m]
        sl[m]={'tasks':len(a),'passed':sum(r['target_pass']=='True' for r in a),
               'items':sum(int(r['items']) for r in a),'undercovered_items':sum(int(r['undercovered_items']) for r in a),
               'energy_evaluations':sum(int(r['energy_evaluations']) for r in a),
               'single_run_seconds':sum(float(r['seconds']) for r in a),
               'maximum_error_to_indicator':max(float(r['maximum_underestimate_ratio']) for r in a)}
    result={'main_24':sm,'adversarial_12_posthoc':ss,'legacy_308_posthoc':sl,
            'claim':'empirical benchmark; no rigorous error bound, no novelty or publication claim'}
    (ROOT/'results/summary.json').write_text(json.dumps(result,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.18,'savefig.dpi':190})
    fig,axs=plt.subplots(1,2,figsize=(11,4.3))
    for ax,rows,title in zip(axs,[main_rows,stress_rows],['Frozen main set (24 cases)','Posthoc adversarial set (12 cases)']):
        for m in METHODS:
            a=[r for r in rows if r['method']==m]
            ax.scatter([r['energy'] for r in a],[max(r['rel']*100,1e-12) for r in a],label=LABELS[m],color=COLORS[m],s=34,alpha=.8)
        ax.axhline(.1,color='k',ls='--',lw=1,label='0.1% target')
        ax.set(xlabel='Real-energy transport evaluations',ylabel='Absolute relative current error (%)',yscale='log',title=title)
        ax.legend(fontsize=7,loc='upper right')
    fig.tight_layout();fig.savefig(ROOT/'results/accuracy_cost.png');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.2));x=np.arange(3)
    axs[0].bar(x,[sl[m]['undercovered_items'] for m in METHODS],color=[COLORS[m] for m in METHODS]);axs[0].set(ylabel='Items with actual error > indicator',title='Legacy replay: all 308 original items',ylim=(0,6))
    axs[1].bar(x,[sl[m]['energy_evaluations']/sl['uniform']['energy_evaluations'] for m in METHODS],color=[COLORS[m] for m in METHODS]);axs[1].axhline(1,color='k',ls='--',lw=1);axs[1].set(ylabel='Transport evaluation count / uniform',title='Identical phase2 goal budgets')
    for ax in axs:ax.set_xticks(x,['Uniform','Spectral','Pole width'])
    for i,m in enumerate(METHODS):
        axs[0].text(i,sl[m]['undercovered_items']+.15,str(sl[m]['undercovered_items']),ha='center')
        axs[1].text(i,sl[m]['energy_evaluations']/sl['uniform']['energy_evaluations']+.04,f"{sl[m]['energy_evaluations']/sl['uniform']['energy_evaluations']:.2f}x",ha='center')
    axs[1].set_ylim(0,2.4);fig.tight_layout();fig.savefig(ROOT/'results/legacy_tradeoff.png');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.1))
    for ax,summary,title in zip(axs,[sm,ss],['Main: 24 cases','Adversarial: 12 cases']):
        base=summary['uniform']['sum_per_case_median_seconds']
        vals=[summary[m]['sum_per_case_median_seconds']/base for m in METHODS]
        ax.bar(x,vals,color=[COLORS[m] for m in METHODS]);ax.axhline(1,color='k',ls='--',lw=1)
        ax.set_xticks(x,['Uniform','Spectral','Pole width']);ax.set(ylabel='Total of case medians / uniform',title=title,ylim=(0,1.3))
        for i,v in enumerate(vals):ax.text(i,v+.03,f'{v:.2f}',ha='center')
    fig.suptitle('Descriptive timings: 3 repetitions, discovery included; no speedup claim',fontsize=11)
    fig.tight_layout();fig.savefig(ROOT/'results/timing.png');plt.close(fig)
    # Narrow physical peak evaluated directly on the real axis.
    case=json.loads((ROOT/'results_stress/cases.json').read_text())[4]
    dev=make_device(case,case['height_eV']);d=discover(dev,case['bias'],case['temperature']);z=min(d['poles'],key=lambda z:z.real);w=-2*z.imag
    offsets=np.linspace(-8,8,1201);energies=z.real+offsets*w;trans=transport(dev,energies,case['bias'])['T']
    fig,ax=plt.subplots(figsize=(8,4));ax.plot(offsets,trans,label='Frozen real-axis transport kernel',color=COLORS['pole'])
    ax.axvline(0,ls=':',color='k',label='Real part of outgoing pole');ax.axvline(-.5,ls='--',color='grey');ax.axvline(.5,ls='--',color='grey')
    ax.set(xlabel='(Energy - pole center) / pole width',ylabel='Transmission probability',title=f"{case['id']}/{case['sample']}: center {z.real:.8f} eV, width {w:.3g} eV")
    ax.legend(fontsize=8);fig.tight_layout();fig.savefig(ROOT/'results/narrow_resonance.png');plt.close(fig)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
