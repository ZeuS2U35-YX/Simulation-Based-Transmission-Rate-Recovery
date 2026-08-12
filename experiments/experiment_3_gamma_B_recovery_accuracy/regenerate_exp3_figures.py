import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=os.path.dirname(__file__)
res=os.path.join(ROOT,'results','combined')
figdir=os.path.join(ROOT,'figures')
os.makedirs(figdir,exist_ok=True)

plt.rcParams.update({
    'font.family':'DejaVu Sans',
    'font.size':9.5,
    'axes.labelsize':11,
    'legend.fontsize':9,
    'xtick.labelsize':9,
    'ytick.labelsize':9,
    'axes.linewidth':0.8,
    'xtick.major.width':0.8,
    'ytick.major.width':0.8,
    'xtick.major.size':4,
    'ytick.major.size':4,
    'pdf.fonttype':42,
    'ps.fonttype':42,
})

BLUE = '#6BAED6'
GRAY_FILL = '0.85'


def polish(ax):
    ax.grid(False)
    ax.tick_params(direction='out', top=False, right=False)
    for s in ax.spines.values():
        s.set_linewidth(0.8)
        s.set_color('black')

paths=pd.read_csv(os.path.join(res,'combined_filtered_B_paths.csv'))
runs=pd.read_csv(os.path.join(res,'combined_mif2_results.csv'))
selected=pd.read_csv(os.path.join(
    ROOT,'results','selected_trajectory','experiment3_task145_B_trajectory.csv'
))

# Recompute summaries from stored combined outputs.
paths['error'] = paths['B_filtered_mean'] - paths['B_true']
paths['sq'] = paths['error']**2
err_rows=[]
for task_id, g in paths.groupby('task_id'):
    err_rows.append({
        'task_id': task_id,
        'RSS': g['sq'].sum(),
        'RMSE': float(np.sqrt(g['sq'].mean())),
        'bias_all': g['error'].mean(),
        'bias_before': g.loc[g['week'] < 5, 'error'].mean(),
        'bias_after': g.loc[g['week'] >= 5, 'error'].mean(),
    })
error_summary = pd.DataFrame(err_rows)

runs = runs[np.isfinite(runs['logLik'])].copy()
gap_rows=[]
for task_id, g in runs.groupby('task_id'):
    z = g.sort_values('logLik', ascending=False).reset_index(drop=True)
    second = z.loc[1, 'logLik'] if len(z) >= 2 else np.nan
    gap_rows.append({
        'task_id': task_id,
        'best_run': int(z.loc[0, 'run']),
        'best_logLik': z.loc[0, 'logLik'],
        'second_best_logLik': second,
        'logLik_gap': z.loc[0, 'logLik'] - second if len(z) >= 2 else np.nan,
    })
likelihood_gaps = pd.DataFrame(gap_rows)

# This optional generator deliberately leaves the canonical metric CSVs unchanged.

# 1. One prespecified task-145 ancestry-preserving particle trajectory.
required={'experiment','task_id','week','B_trajectory','B_true','is_time_zero'}
if not required.issubset(selected.columns) or len(selected) != 71:
    raise ValueError('Invalid selected task-145 trajectory artifact')
if set(selected.task_id) != {145} or selected.week.duplicated().any():
    raise ValueError('Selected trajectory must contain only task 145 at unique times')
fig, ax = plt.subplots(figsize=(6.0, 3.9))
ax.plot(selected['week'], selected['B_trajectory'], color='#0072B2', linewidth=1.45, label='Gamma-noise trajectory (task 145)')
ax.hlines(4.0, xmin=0, xmax=5, color='black', linewidth=1.8, label='True B(t)')
ax.hlines(2.0, xmin=5, xmax=10, color='black', linewidth=1.8)
ax.axvline(5, color='0.72', linewidth=0.8, linestyle=(0, (3, 3)))
ax.set_xlim(0, 10)
ax.set_ylim(min(0,selected.B_trajectory.min())-.2,max(6,selected.B_trajectory.max())+.2)
ax.set_xticks(np.arange(0, 11, 2))
ax.set_yticks(np.arange(0, 7, 1))
ax.set_xlabel('Week')
ax.set_ylabel('Transmission rate, B(t)')
polish(ax)
handles, labels = ax.get_legend_handles_labels(); order = [1, 0]
ax.legend([handles[i] for i in order], [labels[i] for i in order], loc='upper right', frameon=False, handlelength=3.0)
fig.subplots_adjust(left=0.13, right=0.98, bottom=0.17, top=0.98)
fig.savefig(os.path.join(figdir, '01_selected_task_B_trajectory.pdf'))
plt.close(fig)

# 2. RSS distribution
fig, ax = plt.subplots(figsize=(6.0, 3.9))
upper = max(60, int(np.ceil(error_summary['RSS'].max() / 5.0)) * 5)
bins = np.arange(0, upper + 5, 5)
ax.hist(error_summary['RSS'], bins=bins, color=GRAY_FILL, edgecolor='black', linewidth=0.6)
ax.set_xlabel('Observation-time filtering-mean RSS')
ax.set_ylabel('Frequency')
polish(ax)
fig.subplots_adjust(left=0.13, right=0.98, bottom=0.17, top=0.98)
fig.savefig(os.path.join(figdir, 'rss_distribution.pdf'))
plt.close(fig)

# 3. RMSE distribution
fig, ax = plt.subplots(figsize=(6.0, 3.9))
upper = max(1.1, np.ceil(error_summary['RMSE'].max() * 20) / 20)
bins = np.arange(0.3, upper + 0.05, 0.05)
ax.hist(error_summary['RMSE'], bins=bins, color=GRAY_FILL, edgecolor='black', linewidth=0.6)
ax.set_xlabel('Observation-time filtering-mean RMSE')
ax.set_ylabel('Frequency')
polish(ax)
fig.subplots_adjust(left=0.13, right=0.98, bottom=0.17, top=0.98)
fig.savefig(os.path.join(figdir, 'rmse_distribution.pdf'))
plt.close(fig)

# 4. Bias panels
fig, axs = plt.subplots(1, 3, figsize=(7.8, 3.1), sharey=False)
series = [
    ('Overall', error_summary['bias_all'].to_numpy()),
    ('Before week 5', error_summary['bias_before'].to_numpy()),
    ('Week 5 onward', error_summary['bias_after'].to_numpy())
]
allvals = np.concatenate([v for _, v in series])
xmin = min(-1.2, np.floor(allvals.min() / 0.2) * 0.2)
xmax = max(1.2, np.ceil(allvals.max() / 0.2) * 0.2)
bins = np.arange(xmin, xmax + 0.2, 0.2)
for i, (ax, (title, data)) in enumerate(zip(axs, series)):
    ax.hist(data, bins=bins, color=GRAY_FILL, edgecolor='black', linewidth=0.55)
    ax.axvline(0, color='black', linewidth=0.8, linestyle=(0, (4, 3)))
    ax.axvline(np.mean(data), color='0.35', linewidth=0.8, linestyle=(0, (1.2, 2.0)))
    ax.set_title(title, fontsize=9.5, pad=4)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel('Mean estimation error')
    if i == 0:
        ax.set_ylabel('Frequency')
    polish(ax)
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.22, top=0.90, wspace=0.24)
fig.savefig(os.path.join(figdir, 'bias_before_after.pdf'))
plt.close(fig)

# 5. Likelihood-gap distribution
fig, ax = plt.subplots(figsize=(6.0, 3.9))
upper = max(0.15, np.ceil(np.nanmax(likelihood_gaps['logLik_gap']) * 100) / 100)
bins = np.arange(0, upper + 0.005, 0.005)
ax.hist(likelihood_gaps['logLik_gap'].dropna(), bins=bins, color=GRAY_FILL, edgecolor='black', linewidth=0.6)
ax.set_xlabel('Best - second-best evaluated log likelihood')
ax.set_ylabel('Frequency')
polish(ax)
fig.subplots_adjust(left=0.13, right=0.98, bottom=0.22, top=0.98)
fig.savefig(os.path.join(figdir, 'likelihood_gap_distribution.pdf'))
plt.close(fig)

print('Updated figures written to', figdir)
