import os, math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=os.path.dirname(__file__)

plt.rcParams.update({
    'font.family':'DejaVu Sans',
    'font.size':9,
    'axes.labelsize':10,
    'axes.titlesize':9,
    'legend.fontsize':8.5,
    'xtick.labelsize':8.5,
    'ytick.labelsize':8.5,
    'axes.linewidth':0.8,
    'xtick.major.width':0.8,
    'ytick.major.width':0.8,
    'xtick.major.size':4,
    'ytick.major.size':4,
    'pdf.fonttype':42,
    'ps.fonttype':42,
})

def polish(ax):
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.8); s.set_color('black')
    ax.tick_params(direction='out', top=False, right=False)
    ax.grid(False)

def get_data(combined_root, comparison_root):
    paired=pd.read_csv(os.path.join(comparison_root,'paired_model_comparison.csv'))
    gp=pd.read_csv(os.path.join(combined_root,'gamma','combined_B_paths.csv'))
    cp=pd.read_csv(os.path.join(combined_root,'constant','combined_B_paths.csv'))
    return paired,gp,cp

def mean_paths(gp,cp):
    mg=gp.groupby('week',as_index=False)[['B_estimate','B_true']].mean()
    mc=cp.groupby('week',as_index=False)['B_estimate'].mean()
    return mg,mc

def save_mean_path(gp,cp,out):
    mg,mc=mean_paths(gp,cp)
    # restrained publication colors: true path in solid black, gamma in light blue,
    # constant-B in brown-red; gray switch-time reference retained.
    true_col='black'
    gamma_col='#5DA5DA'
    const_col='#A05A4A'
    switch_col='0.70'
    fig,ax=plt.subplots(figsize=(6.15,4.15))
    # Draw the true B(t) explicitly as a piecewise-constant function rather than
    # as one connected curve.
    ax.hlines(4.0, xmin=0, xmax=5, color=true_col, linewidth=1.6, label='True B(t)')
    ax.hlines(2.0, xmin=5, xmax=10, color=true_col, linewidth=1.6)
    ax.plot(mg.week, mg.B_estimate, color=gamma_col, linewidth=1.5, label='Gamma-noise model')
    ax.plot(mc.week, mc.B_estimate, color=const_col, linewidth=1.35, linestyle=(0,(1.5,2.2)), label='Constant-B')
    ax.axvline(5, color=switch_col, linewidth=0.8, linestyle=(0,(3,3)))
    ax.set_xlim(0,10); ax.set_ylim(0,6)
    ax.set_xticks(np.arange(0,11,2)); ax.set_yticks(np.arange(0,7,1))
    ax.set_xlabel('Week'); ax.set_ylabel('Transmission rate, B(t)')
    polish(ax)
    ax.legend(loc='upper right', frameon=False, handlelength=3.2)
    fig.subplots_adjust(left=.13,right=.98,bottom=.15,top=.98)
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)

def stacked_hist(a,b,out,xlabel,labels=('Gamma-noise model','Constant-B'),bins=None,xlim=None,vline=None):
    if bins is None:
        lo=min(np.min(a),np.min(b)); hi=max(np.max(a),np.max(b)); bins=np.linspace(lo,hi,16)
    fig,axs=plt.subplots(2,1,figsize=(6.15,4.6),sharex=True)
    for ax,data,label in zip(axs,[a,b],labels):
        ax.hist(data,bins=bins,color='0.82',edgecolor='black',linewidth=0.75)
        if vline is not None: ax.axvline(vline,color='black',linewidth=0.8,linestyle=(0,(4,3)))
        polish(ax)
        ax.text(.98,.86,label,transform=ax.transAxes,ha='right',va='top',fontsize=9)
        ax.set_ylabel('Frequency')
    axs[1].set_xlabel(xlabel)
    if xlim: axs[1].set_xlim(*xlim)
    fig.subplots_adjust(left=.13,right=.98,bottom=.13,top=.98,hspace=.10)
    fig.savefig(out,bbox_inches='tight')
    plt.close(fig)

def save_bias(p,out):
    series=[
        ('Overall',p.bias_all_gamma,p.bias_all_constant),
        ('Before week 5',p.bias_before_gamma,p.bias_before_constant),
        ('From week 5',p.bias_after_gamma,p.bias_after_constant),
    ]
    xmin=-1.25; xmax=2.75; bins=np.arange(xmin,xmax+0.2001,0.2)
    fig,axs=plt.subplots(2,3,figsize=(7.15,4.55),sharex=True)
    for j,(title,g,c) in enumerate(series):
        for i,(data,rowlab) in enumerate([(g,'Gamma-noise model'),(c,'Constant-B')]):
            ax=axs[i,j]
            ax.hist(data,bins=bins,color='0.82',edgecolor='black',linewidth=0.65)
            ax.axvline(0,color='black',linewidth=.8,linestyle=(0,(4,3)))
            ax.axvline(float(np.mean(data)),color='0.35',linewidth=.9,linestyle=(0,(1.2,2.0)))
            polish(ax)
            if i==0: ax.set_title(title,pad=5)
            if j==0: ax.set_ylabel('Frequency')
            if i==1: ax.set_xlabel('Mean estimation error')
            if j==2: ax.text(1.04,.5,rowlab,transform=ax.transAxes,rotation=-90,ha='left',va='center',fontsize=9)
            ax.set_xlim(xmin,xmax)
    fig.subplots_adjust(left=.09,right=.94,bottom=.15,top=.94,wspace=.16,hspace=.18)
    fig.savefig(out,bbox_inches='tight')
    plt.close(fig)

def save_rmse_scatter(p,out):
    fig,ax=plt.subplots(figsize=(4.85,4.75))
    lo=0.30; hi=2.00
    ax.scatter(p.RMSE_constant,p.RMSE_gamma,s=25,facecolors='none',edgecolors='black',linewidths=.8)
    ax.plot([lo,hi],[lo,hi],color='black',linewidth=.9,linestyle=(0,(5,3)))
    ax.set_xlim(lo,hi); ax.set_ylim(lo,hi); ax.set_aspect('equal',adjustable='box')
    ax.set_xlabel('Constant-B RMSE'); ax.set_ylabel('Gamma-noise model RMSE')
    polish(ax)
    fig.subplots_adjust(left=.17,right=.98,bottom=.15,top=.98)
    fig.savefig(out,bbox_inches='tight')
    plt.close(fig)

def save_loglik(p,out):
    x=p.delta_logLik_gamma_minus_constant
    bins=np.arange(-5,46,3)
    fig,ax=plt.subplots(figsize=(6.15,4.0))
    ax.hist(x,bins=bins,color='0.82',edgecolor='black',linewidth=.75)
    ax.axvline(0,color='black',linewidth=.9,linestyle=(0,(5,3)))
    ax.set_xlabel('Independent log likelihood difference (Gamma-noise - constant-B)')
    ax.set_ylabel('Frequency')
    polish(ax)
    fig.subplots_adjust(left=.13,right=.98,bottom=.17,top=.98)
    fig.savefig(out,bbox_inches='tight')
    plt.close(fig)

def generate_compare(combined_root,comparison_root,figdir):
    os.makedirs(figdir,exist_ok=True)
    p,gp,cp=get_data(combined_root,comparison_root)
    save_mean_path(gp,cp,os.path.join(figdir,'01_mean_recovered_B_paths.pdf'))
    bins=np.arange(0,261,15)
    stacked_hist(p.RSS_gamma,p.RSS_constant,os.path.join(figdir,'02_RSS_distributions.pdf'),'Residual sum of squares of B(t)',bins=bins,xlim=(0,260))
    save_bias(p,os.path.join(figdir,'03_bias_distributions.pdf'))
    save_rmse_scatter(p,os.path.join(figdir,'04_paired_RMSE_scatter.pdf'))
    bins=np.arange(.3,2.01,.1)
    stacked_hist(p.RMSE_gamma,p.RMSE_constant,os.path.join(figdir,'05_RMSE_distributions.pdf'),'RMSE of B(t)',bins=bins,xlim=(.3,2.0))
    save_loglik(p,os.path.join(figdir,'06_independent_loglik_difference.pdf'))

def best_runs(runs):
    return runs.loc[runs.groupby('task_id')['logLik'].idxmax(),['task_id','run']].set_index('task_id')['run'].to_dict()

def convergence_pdf(trace_file,runs_file,out,model):
    tr=pd.read_csv(trace_file); runs=pd.read_csv(runs_file); best=best_runs(runs)
    with PdfPages(out) as pdf:
        for task in sorted(tr.task_id.unique()):
            x=tr[tr.task_id==task]
            if model=='gamma': vars=[('loglik','Internal log likelihood'),('B0','B0'),('sigma_beta','sigma_beta')]
            else: vars=[('loglik','Internal log likelihood'),('Beta','Beta')]
            fig,axs=plt.subplots(1,len(vars),figsize=(7.2,2.85 if model=='gamma' else 3.05))
            if len(vars)==1: axs=[axs]
            br=best.get(task,None)
            for ax,(col,lab) in zip(axs,vars):
                for run in sorted(x.run.unique()):
                    z=x[x.run==run].sort_values('mif_iteration')
                    if run==br:
                        ax.plot(z.mif_iteration,z[col],color='black',linewidth=1.15)
                    else:
                        ax.plot(z.mif_iteration,z[col],color='0.68',linewidth=.55)
                ax.axvline(500,color='0.35',linewidth=.75,linestyle=(0,(3,3)))
                ax.set_xlim(0,600); ax.set_xlabel('MIF iteration'); ax.set_ylabel(lab); polish(ax)
            axs[0].text(.02,.96,f'Task {task}',transform=axs[0].transAxes,ha='left',va='top',fontsize=9)
            axs[-1].plot([],[],color='black',linewidth=1.15,label='Best run')
            axs[-1].plot([],[],color='0.68',linewidth=.55,label='Other starts')
            axs[-1].legend(frameon=False,loc='best',fontsize=7.5)
            fig.subplots_adjust(left=.08,right=.99,bottom=.22,top=.96,wspace=.33)
            pdf.savefig(fig,bbox_inches='tight'); plt.close(fig)

final_comb=os.path.join(ROOT,'results','combined')
final_comp=os.path.join(ROOT,'results','comparison')
generate_compare(final_comb,final_comp,os.path.join(ROOT,'figures','comparison'))
convergence_pdf(os.path.join(final_comb,'gamma','combined_mif2_traces.csv'),os.path.join(final_comb,'gamma','combined_mif2_results.csv'),os.path.join(ROOT,'figures','convergence','gamma_convergence_diagnostics.pdf'),'gamma')
convergence_pdf(os.path.join(final_comb,'constant','combined_mif2_traces.csv'),os.path.join(final_comb,'constant','combined_mif2_results.csv'),os.path.join(ROOT,'figures','convergence','constant_B_convergence_diagnostics.pdf'),'constant')

pilot_comb=os.path.join(ROOT,'results','pilot_combined')
pilot_comp=os.path.join(ROOT,'results','pilot_comparison')
# Pilot outputs are retained unchanged by default. Set this explicit opt-in only
# when the pilot figures themselves are intentionally being regenerated.
if os.environ.get('REGENERATE_EXP4_PILOT_FIGURES') == '1' and os.path.exists(pilot_comp):
    generate_compare(pilot_comb,pilot_comp,os.path.join(ROOT,'figures','pilot'))
    convergence_pdf(os.path.join(pilot_comb,'gamma','combined_mif2_traces.csv'),os.path.join(pilot_comb,'gamma','combined_mif2_results.csv'),os.path.join(ROOT,'figures','pilot','gamma_convergence_diagnostics.pdf'),'gamma')
    convergence_pdf(os.path.join(pilot_comb,'constant','combined_mif2_traces.csv'),os.path.join(pilot_comb,'constant','combined_mif2_results.csv'),os.path.join(ROOT,'figures','pilot','constant_B_convergence_diagnostics.pdf'),'constant')
print('done')
