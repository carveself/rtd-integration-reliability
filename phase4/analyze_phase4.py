from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
METHODS=['uniform','spectral','pole','selective']
COLORS=['#777777','#0072B2','#D55E00','#009E73']

def read(p):return list(csv.DictReader(p.open()))

def main():
    fresh=read(P/'results/benchmark.csv');legacy=read(P/'results_regression/benchmark.csv')
    summary={'fresh':{},'legacy':{}}
    for m in METHODS:
        a=[r for r in fresh if r['method']==m];b=[r for r in legacy if r['method']==m]
        summary['fresh'][m]={'cases':len(a),'accepted_reference':sum(r['reference_accepted']=='True' for r in a),
           'passed':sum(r['passed']=='True' for r in a),'undercovered':sum(r['covered']=='False' for r in a),
           'maximum_relative_error':max(float(r['relative_error']) for r in a),
           'energy_evaluations':sum(int(r['energy_evaluations']) for r in a),
           'sum_case_median_seconds':sum(float(r['seconds_median']) for r in a)}
        summary['legacy'][m]={'tasks':len(b),'passed':sum(r['passed']=='True' for r in b),
           'undercovered':sum(int(r['undercovered_items']) for r in b),'energy_evaluations':sum(int(r['energy_evaluations']) for r in b),
           'single_run_seconds':sum(float(r['seconds']) for r in b)}
    selective=[r for r in fresh if r['method']=='selective']
    summary['selection']={'eigen_seeds':sum(int(r['eigen_seeds']) for r in selective),
                          'retained_seeds':sum(int(r['retained_seeds']) for r in selective),
                          'coarse_evaluations':sum(int(r['coarse_evaluations']) for r in selective)}
    for name in ['fresh','legacy']:
        uniform=summary[name]['uniform']['energy_evaluations']
        spectral=summary[name]['spectral']['energy_evaluations']
        for m in METHODS:
            s=summary[name][m];s['energy_ratio_to_uniform']=s['energy_evaluations']/uniform
            s['energy_ratio_to_spectral']=s['energy_evaluations']/spectral
    errors=read(P/'results_regression/sample_errors.csv')
    failed=[r for r in errors if r['method']=='selective' and float(r['actual_error_A_m2'])>float(r['estimated_error_A_m2'])]
    for r in failed:r['actual_to_estimate_ratio']=float(r['actual_error_A_m2'])/float(r['estimated_error_A_m2'])
    summary['selective_regression_failures']=failed
    summary['verdict']='Selective method saves evaluations but does not preserve estimator coverage on old cases. Experimental, not a reliable default or publication-ready new method.'
    (P/'results/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.18,'savefig.dpi':180})
    fig,axs=plt.subplots(1,2,figsize=(11,4.2));x=np.arange(4)
    axs[0].bar(x,[summary['fresh'][m]['energy_ratio_to_uniform'] for m in METHODS],color=COLORS)
    axs[0].axhline(1,ls='--',color='k');axs[0].set(ylabel='All real-energy evaluations / uniform',title='Fresh 24 cases: all meet 0.1% current target',ylim=(0,1.8))
    axs[1].bar(x,[summary['legacy'][m]['undercovered'] for m in METHODS],color=COLORS)
    axs[1].set(ylabel='Actual error exceeds empirical indicator',title='Legacy 308 items: same goal budgets',ylim=(0,6))
    for ax in axs:ax.set_xticks(x,['Uniform','Spectral','Pole','Selective'])
    for i,m in enumerate(METHODS):
        axs[0].text(i,summary['fresh'][m]['energy_ratio_to_uniform']+.035,f"{summary['fresh'][m]['energy_ratio_to_uniform']:.2f}x",ha='center')
        axs[1].text(i,summary['legacy'][m]['undercovered']+.12,str(summary['legacy'][m]['undercovered']),ha='center')
    fig.tight_layout();fig.savefig(P/'results/efficiency_reliability.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.4))
    for m,color in zip(METHODS,COLORS):
        a=[r for r in fresh if r['method']==m]
        ax.scatter([int(r['energy_evaluations']) for r in a],[max(float(r['relative_error'])*100,1e-12) for r in a],label=m,color=color,s=35,alpha=.75)
    ax.axhline(.1,ls='--',color='k',label='0.1% target');ax.set(xlabel='Energy evaluations (coarse probe included)',ylabel='Current relative error (%)',yscale='log',title='Fresh validation: accuracy alone does not test old coverage failures');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(P/'results/fresh_error_cost.png');plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
