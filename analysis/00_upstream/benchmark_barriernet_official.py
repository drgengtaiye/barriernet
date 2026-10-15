#!/usr/bin/env python3
"""
=======================================================================
PRJEB25780 ICI Response Benchmark (using barriernet official genesets)
=======================================================================
用法:
    pip install tidepy scikit-learn matplotlib pandas numpy scipy pyyaml
    python benchmark_barriernet_official.py

输入: mGC_Pembrolizumab_Kim2018.TPM (和本脚本同目录)
      geneset_rna.yaml (可选, 不提供则用内嵌定义)

输出: Fig_benchmark_main.pdf/png   — 主图 (waterfall + boxplot + ROC)
      Fig_AUC_bar.pdf/png          — AUC对比柱状图
      Fig_heatmap.pdf/png          — 患者×score热图
      Fig_MSS_subgroup.pdf/png     — MSS亚组
      benchmark_scores.tsv         — 完整score表
=======================================================================
"""
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
matplotlib.rcParams.update({
    'pdf.fonttype': 42, 'font.family': 'Arial', 'font.size': 10,
    'axes.linewidth': 0.8,
})
import warnings; warnings.filterwarnings('ignore')
import os, sys

# ======================= CONFIG =======================
TPM_FILE = 'mGC_Pembrolizumab_Kim2018.TPM'
YAML_FILE = 'geneset_rna.yaml'  # optional
# ======================================================

# =================================================================
# PART A: barriernet geneset (hardcoded from geneset_rna.yaml v1)
#         So this script works WITHOUT installing barriernet
# =================================================================
BARRIERNET_MODULES = {
    'CAF_ECM': ['COL1A1','COL1A2','LUM','DCN','POSTN','TAGLN','ACTA2','THY1','FAP','FN1','COL3A1'],
    'TAM_SUP': ['SPP1','APOE','TREM2','C1QA','C1QB','C1QC','LGALS3','MARCO','IL10','SOCS3','STAT3'],
    'APM':     ['B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2','PSMB8','PSMB9'],
    'IFN':     ['STAT1','IRF1','CXCL9','CXCL10','IFNG'],
    'TLS_B':   ['CXCL13','LTB','MS4A1','CD79A','CD74','MZB1'],
    'T_CELL':  ['CD3D','CD3E','TRAC','TRBC1','CD8A','GZMB','NKG7'],
}
BARRIERNET_AXES = {
    'Barrier':   {'positive_modules': ['CAF_ECM','TAM_SUP']},
    'Immune':    {'positive_modules': ['APM','IFN','TLS_B','T_CELL']},
    'NetImmune': {'positive_modules': ['APM','IFN','TLS_B','T_CELL'],
                  'negative_modules': ['CAF_ECM','TAM_SUP']},
}

# Try loading from YAML if available (overrides hardcoded)
try:
    import yaml
    if os.path.exists(YAML_FILE):
        with open(YAML_FILE) as f:
            obj = yaml.safe_load(f)
        BARRIERNET_MODULES = {k: [str(x).upper() for x in v] for k,v in obj['modules'].items()}
        BARRIERNET_AXES = obj['axes']
        print(f"✓ Loaded geneset from {YAML_FILE}")
    else:
        print(f"  {YAML_FILE} not found, using built-in barriernet geneset_rna v1")
except ImportError:
    print("  pyyaml not installed, using built-in barriernet geneset_rna v1")

# =================================================================
# PART B: barriernet core scoring (exact replica of core.py)
# =================================================================
def zscore_series(x):
    x = pd.to_numeric(x, errors="coerce")
    mu = np.nanmean(x.values)
    sd = np.nanstd(x.values, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * 0.0
    return (x - mu) / (sd + 1e-9)

def score_modules(expr, modules, min_genes=3):
    """Exact barriernet core.score_modules logic."""
    X = expr.copy()
    X.columns = [str(c).upper() for c in X.columns]
    X = X.apply(zscore_series, axis=0)  # z-score per gene across samples

    scores = pd.DataFrame(index=X.index)
    used = {}
    for m, genes in modules.items():
        g = [str(x).upper() for x in genes]
        g2 = [x for x in g if x in X.columns]
        used[m] = g2
        if len(g2) < min_genes:
            scores[m] = np.nan
            print(f"  ⚠ {m}: only {len(g2)}/{len(g)} genes found (< min_genes={min_genes})")
        else:
            scores[m] = X[g2].mean(axis=1, skipna=True)
            print(f"  ✓ {m}: {len(g2)}/{len(g)} genes")
    return scores, used

def score_axes(module_scores, axes):
    """Exact barriernet core.score_axes logic."""
    out = pd.DataFrame(index=module_scores.index)
    for ax, spec in axes.items():
        pos = [x for x in spec.get('positive_modules', []) if x in module_scores.columns]
        neg = [x for x in spec.get('negative_modules', []) if x in module_scores.columns]
        pos_v = module_scores[pos].mean(axis=1, skipna=True) if pos else 0.0
        neg_v = module_scores[neg].mean(axis=1, skipna=True) if neg else 0.0
        out[ax] = pos_v - neg_v
    return out

# =================================================================
# PART C: Clinical data (Kim 2018 Table 2, all 61 patients)
# =================================================================
CLINICAL_RAW = """EP-36,PB-16-045,MSI,Negative,CR
EP-60,PB-16-068,MSI,Negative,CR
EP-32,PB-16-031,MSS,Negative,CR
EP-05,PB-16-005,MSS,Negative,PR
EP-10,PB-16-010,MSI,Negative,PR
EP-11,PB-16-013,MSS,Negative,PR
EP-23,PB-16-019,MSI,Negative,PR
EP-19,PB-16-020,MSS,Positive,PR
EP-24,PB-16-021,MSS,Positive,PR
EP-39,PB-16-038,MSS,Positive,PR
EP-40,PB-16-044,MSI,Negative,PR
EP-43,PB-16-048,MSI,Negative,PR
EP-49,PB-16-053,MSS,Positive,PR
EP-54,PB-16-059,MSS,Positive,PR
EP-56,PB-16-063,MSS,Positive,PR
EP-58,PB-16-066,MSS,Negative,SD
EP-06,PB-16-006,MSS,Negative,SD
EP-13,PB-16-016,MSS,Negative,SD
EP-18,PB-16-022,MSS,Negative,SD
EP-27,PB-16-026,MSS,Negative,SD
EP-29,PB-16-034,MSS,Negative,SD
EP-31,PB-16-037,MSS,Negative,SD
EP-34,PB-16-040,MSS,Negative,SD
EP-35,PB-16-042,MSS,Negative,SD
EP-42,PB-16-047,MSS,Negative,SD
EP-41,PB-16-049,MSS,Negative,SD
EP-50,PB-16-052,MSS,Negative,SD
EP-47,PB-16-054,MSS,Negative,SD
EP-45,PB-16-056,MSS,Negative,SD
EP-52,PB-16-058,MSS,Negative,SD
EP-57,PB-16-062,MSS,Negative,SD
EP-55,PB-16-064,MSS,Negative,SD
EP-61,PB-16-069,MSS,Negative,SD
EP-03,PB-16-003,MSS,Negative,SD
EP-12,PB-16-014,MSS,Negative,SD
EP-01,PB-16-001,MSS,Negative,PD
EP-02,PB-16-002,MSS,Negative,PD
EP-04,PB-16-004,MSS,Negative,PD
EP-07,PB-16-011,MSS,Negative,PD
EP-14,PB-16-015,MSS,Negative,PD
EP-16,PB-16-018,MSS,Negative,PD
EP-17,PB-16-023,MSS,Negative,PD
EP-21,PB-16-024,MSS,Negative,PD
EP-20,PB-16-025,MSS,Negative,PD
EP-28,PB-16-029,MSS,Negative,PD
EP-22,PB-16-030,MSS,Negative,PD
EP-25,PB-16-032,MSS,Negative,PD
EP-30,PB-16-035,MSS,Negative,PD
EP-33,PB-16-039,MSS,Negative,PD
EP-37,PB-16-041,MSS,Negative,PD
EP-38,PB-16-043,MSS,Negative,PD
EP-44,PB-16-051,MSS,Negative,PD
EP-48,PB-16-055,MSI,Negative,PD
EP-46,PB-16-057,MSS,Negative,PD
EP-53,PB-16-060,MSS,Negative,PD
EP-59,PB-16-067,MSS,Negative,PD
EP-51,PB-19-009,MSS,Negative,PD
EP-09,PB-16-007,MSS,Negative,NE
EP-08,PB-16-008,MSS,Negative,NE
EP-15,PB-16-017,MSS,Negative,NE
EP-26,PB-16-033,MSS,Negative,NE"""

import io
clin = pd.read_csv(io.StringIO(CLINICAL_RAW), header=None,
    names=['Patient_ID','SEQ_ID','MSI_pentaplex','EBV_insitu','Response'])
clin['MSI_H'] = clin['MSI_pentaplex'] == 'MSI'
clin['EBV_pos'] = clin['EBV_insitu'] == 'Positive'
clin['Responder'] = clin['Response'].isin(['CR','PR']).astype(int)
clin['Evaluable'] = clin['Response'] != 'NE'

# =================================================================
# PART D: Load TPM + merge + score
# =================================================================
print("=" * 65)
print("  PRJEB25780 Benchmark (barriernet official geneset)")
print("=" * 65)

if not os.path.exists(TPM_FILE):
    sys.exit(f"ERROR: {TPM_FILE} not found")

with open(TPM_FILE) as f:
    header = f.readline().strip().split('\t')
expr = pd.read_csv(TPM_FILE, sep='\t', header=None, skiprows=1, index_col=0)
expr.columns = header
expr = expr.T.astype(float)
print(f"\n✓ TPM loaded: {expr.shape[0]} samples × {expr.shape[1]} genes")

# Merge with clinical
clin_idx = clin.set_index('SEQ_ID')
merged_clin = clin_idx.loc[clin_idx.index.isin(expr.index)]
expr_matched = expr.loc[merged_clin.index]
evaluable = merged_clin['Evaluable']
expr_eval = expr_matched[evaluable]
clin_eval = merged_clin[evaluable]
y = clin_eval['Responder'].values
print(f"✓ Matched: {len(expr_eval)} evaluable ({y.sum()} R, {(1-y).sum()} NR)")

# --- barriernet scores ---
print(f"\n--- barriernet modules ---")
mod_scores, used = score_modules(expr_eval, BARRIERNET_MODULES, min_genes=3)
axes_scores = score_axes(mod_scores, BARRIERNET_AXES)

print(f"\n--- barriernet axes ---")
for ax in axes_scores.columns:
    print(f"  {ax}: mean={axes_scores[ax].mean():.3f}, std={axes_scores[ax].std():.3f}")

# Z-standardize axes
for c in axes_scores.columns:
    axes_scores[f'{c}_z'] = zscore_series(axes_scores[c])

# =================================================================
# PART E: Benchmark scores (GEP-18, IPS, IFNG, TIDE)
# =================================================================
print(f"\n--- Benchmark methods ---")
X_upper = expr_eval.copy()
X_upper.columns = [str(c).upper() for c in X_upper.columns]
X_z = X_upper.apply(zscore_series, axis=0)

def zmod(genes, name=''):
    avail = [g for g in genes if g in X_z.columns]
    if name:
        miss = [g for g in genes if g not in X_z.columns]
        tag = f"({len(avail)}/{len(genes)})" + (f" miss:{miss}" if miss else " ✓")
        print(f"  {name}: {tag}")
    return X_z[avail].mean(axis=1) if len(avail) >= 2 else pd.Series(0, index=X_z.index)

# GEP-18 (Cristescu/Merck pembrolizumab biomarker)
GEP18 = zmod(['CCL5','CD27','CD274','CD276','CD8A','CMKLR1','CXCL9','CXCR6',
    'HLA-DQA1','HLA-DRB1','HLA-E','IDO1','LAG3','NKG7','PDCD1LG2','PSMB10','STAT1','TIGIT'], 'GEP-18')

# IPS (Charoentong 2017)
ips_mhc = zmod(['HLA-A','HLA-B','HLA-C','HLA-DPA1','HLA-DPB1','HLA-DRA','TAP1','TAP2'], 'IPS_MHC')
ips_ec  = zmod(['CD8A','CD8B','GZMA','GZMB','PRF1','IFNG'], 'IPS_EC')
ips_sc  = zmod(['CTLA4','PDCD1','LAG3','TIGIT','HAVCR2'], 'IPS_SC')
ips_cp  = zmod(['CD274','PDCD1LG2'], 'IPS_CP')
IPS = ips_mhc + ips_ec - ips_sc - ips_cp

# IFNG 6-gene
IFNG_score = zmod(['IFNG','CXCL9','CXCL10','IDO1','STAT1','HLA-DRA'], 'IFNG_6gene')

# PD-L1
PDL1 = X_upper['CD274'].values if 'CD274' in X_upper.columns else np.zeros(len(expr_eval))

# TIDE
TIDE_score = None
try:
    from tidepy.pred import TIDE as TIDEpred
    expr_for_tide = expr_eval.T
    expr_norm = expr_for_tide.subtract(expr_for_tide.mean(axis=1), axis=0)
    tide_result = TIDEpred(expr_norm, cancer='Other', pretreat=True, vthres=0.)
    tide_result = tide_result.loc[expr_eval.index]
    TIDE_score = tide_result['TIDE'].values
    print(f"  TIDE: ✓ (via tidepy)")
except ImportError:
    print(f"  TIDE: tidepy not installed, using exclusion proxy")
except Exception as e:
    print(f"  TIDE: error ({e}), using exclusion proxy")

if TIDE_score is None:
    tide_excl = (zmod(['COL1A1','COL1A2','COL3A1','FN1','POSTN','FAP','ACTA2','PDGFRB']) +
                 zmod(['CD163','MSR1','MRC1','IL10','TGFB1','SPP1','APOE','C1QA'])) / 2
    TIDE_score = tide_excl.values

# Clinical
MSI_H = clin_eval['MSI_H'].astype(int).values
EBV_pos = clin_eval['EBV_pos'].astype(int).values
MSI_EBV = ((clin_eval['MSI_H']) | (clin_eval['EBV_pos'])).astype(int).values

# =================================================================
# PART F: Build score table + compute AUCs
# =================================================================
scores = pd.DataFrame({
    'Patient_ID': clin_eval['Patient_ID'].values,
    'Response': clin_eval['Response'].values,
    'Responder': y,
    'MSI_H': clin_eval['MSI_H'].values,
    'EBV_pos': clin_eval['EBV_pos'].values,
    # barriernet modules
    'CAF_ECM': mod_scores['CAF_ECM'].values,
    'TAM_SUP': mod_scores['TAM_SUP'].values,
    'APM': mod_scores['APM'].values,
    'IFN': mod_scores['IFN'].values,
    'TLS_B': mod_scores['TLS_B'].values,
    'T_CELL': mod_scores['T_CELL'].values,
    # barriernet axes
    'Barrier': axes_scores['Barrier'].values,
    'Immune': axes_scores['Immune'].values,
    'NetImmune': axes_scores['NetImmune'].values,
    'NetImmune_z': axes_scores['NetImmune_z'].values,
    'Barrier_z': axes_scores['Barrier_z'].values,
    # benchmarks
    'GEP18': GEP18.values,
    'IPS': IPS.values,
    'IFNG_score': IFNG_score.values,
    'PDL1_expr': PDL1,
    'TIDE': TIDE_score,
}, index=expr_eval.index)

# Combined
scores['NetImm_plus_MSI_EBV'] = scores['NetImmune_z'].values + 3 * MSI_EBV

scores.to_csv('benchmark_scores.tsv', sep='\t')
print(f"\n✓ Saved: benchmark_scores.tsv")

# --- AUC table ---
predictors = {
    'NetImmune':      scores['NetImmune_z'].values,
    'Barrier(neg)':  -scores['Barrier_z'].values,
    'GEP-18':         scores['GEP18'].values,
    'IPS':            scores['IPS'].values,
    'IFNG':           scores['IFNG_score'].values,
    'PD-L1':          scores['PDL1_expr'],
    'TIDE(neg)':     -scores['TIDE'].values,
    'MSI-H':          MSI_H.astype(float),
    'EBV+':           EBV_pos.astype(float),
    'MSI-H|EBV+':     MSI_EBV.astype(float),
    'NetImm+MSI|EBV': scores['NetImm_plus_MSI_EBV'].values,
}

print(f"\n{'='*58}")
print(f"  {'Method':<20} {'AUC':>6} {'p(MWU)':>10} {'Cliff_d':>8}")
print(f"  {'-'*52}")
auc_dict = {}
for name, ys in predictors.items():
    auc = roc_auc_score(y, ys)
    r_s, nr_s = ys[y==1], ys[y==0]
    try: _, p = stats.mannwhitneyu(r_s, nr_s, alternative='two-sided')
    except: p = 1.0
    n1, n2 = len(r_s), len(nr_s)
    cd = (sum(1 for a in r_s for b in nr_s if a>b) -
          sum(1 for a in r_s for b in nr_s if a<b)) / (n1*n2) if n1*n2 else 0
    auc_dict[name] = auc
    sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
    ours = ' ← OURS' if 'Net' in name or 'Barrier' in name else ''
    print(f"  {name:<20} {auc:>6.3f} {p:>10.4f} {cd:>+8.3f}  {sig}{ours}")
print(f"  {'='*52}")

# =================================================================
# PART G: FIGURES
# =================================================================
print("\n--- Generating figures ---")
C = {'red':'#E74C3C','blue':'#3498DB','purple':'#9B59B6','orange':'#F39C12',
     'green':'#2ECC71','gray':'#95A5A6','teal':'#1ABC9C','dark':'#2C3E50'}
resp_colors = {'CR':'#1a9641','PR':'#a6d96a','SD':'#fdae61','PD':'#d7191c'}

# ============ FIG 1: Main 3-panel ============
fig = plt.figure(figsize=(16, 5))
gs = gridspec.GridSpec(1, 3, wspace=0.3, width_ratios=[1.2, 1, 1])

# --- A: Boxplots ---
ax = fig.add_subplot(gs[0])
box_methods = [('NetImmune','NetImmune_z'),('GEP-18','GEP18'),('IFNG','IFNG_score'),
               ('TIDE(neg)','TIDE'),('IPS','IPS')]
for i, (label, col) in enumerate(box_methods):
    vals = -scores[col].values if 'TIDE' in label else scores[col].values
    r_v, nr_v = vals[y==1], vals[y==0]
    _, p = stats.mannwhitneyu(r_v, nr_v, alternative='two-sided')
    pos_nr, pos_r = i*2.8, i*2.8+1
    bp = ax.boxplot([nr_v, r_v], positions=[pos_nr, pos_r], widths=0.6,
                    patch_artist=True, medianprops=dict(color='black', lw=1.5),
                    flierprops=dict(markersize=3))
    bp['boxes'][0].set_facecolor('#D5D8DC'); bp['boxes'][1].set_facecolor('#5DADE2')
    for j, d in enumerate([nr_v, r_v]):
        x = np.random.normal([pos_nr, pos_r][j], 0.05, len(d))
        ax.scatter(x, d, c='black', s=12, alpha=0.4, zorder=5)
    sig = '***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'
    ymax = max(np.max(r_v), np.max(nr_v))
    ax.text((pos_nr+pos_r)/2, ymax+0.15, sig, ha='center', fontsize=9, fontweight='bold')
ax.set_xticks([i*2.8+0.5 for i in range(len(box_methods))])
ax.set_xticklabels([m[0] for m in box_methods], fontsize=9)
ax.set_ylabel('Score')
ax.set_title('A  Responder vs Non-responder', fontsize=11, fontweight='bold', loc='left')

# --- B: ROC ---
ax = fig.add_subplot(gs[1])
roc_items = [
    ('NetImmune', scores['NetImmune_z'].values, C['red'], '-', 2.5),
    ('GEP-18', scores['GEP18'].values, C['purple'], '-', 2),
    ('IFNG', scores['IFNG_score'].values, C['orange'], '-', 2),
    ('TIDE(neg)', -scores['TIDE'].values, C['dark'], '-', 2),
    ('IPS', scores['IPS'].values, C['green'], '-', 2),
    ('PD-L1', PDL1, C['gray'], '-', 1.5),
    ('MSI|EBV', MSI_EBV.astype(float), C['teal'], '--', 1.5),
    ('NetImm+MSI|EBV', scores['NetImm_plus_MSI_EBV'].values, C['red'], '--', 2.5),
]
for name, vals, color, ls, lw in roc_items:
    auc = roc_auc_score(y, vals)
    fpr, tpr, _ = roc_curve(y, vals)
    ax.plot(fpr, tpr, color=color, ls=ls, lw=lw, label=f'{name} ({auc:.3f})')
ax.plot([0,1],[0,1],'k--',lw=0.6,alpha=0.3)
ax.set_xlabel('FPR'); ax.set_ylabel('TPR')
ax.legend(loc='lower right', fontsize=7, framealpha=0.9)
ax.set_title('B  ROC: ICI Response', fontsize=11, fontweight='bold', loc='left')

# --- C: AUC bar ---
ax = fig.add_subplot(gs[2])
bar_items = {k:v for k,v in auc_dict.items() if k not in ['EBV+','NetImm+MSI|EBV']}
sorted_it = sorted(bar_items.items(), key=lambda x: x[1])
names_bar = [x[0] for x in sorted_it]; aucs_bar = [x[1] for x in sorted_it]
colors_bar = []
for n in names_bar:
    if 'Net' in n: colors_bar.append(C['red'])
    elif 'MSI' in n or 'EBV' in n: colors_bar.append(C['teal'])
    elif 'GEP' in n: colors_bar.append(C['purple'])
    elif 'IFNG' in n: colors_bar.append(C['orange'])
    elif 'TIDE' in n: colors_bar.append(C['dark'])
    elif 'Barrier' in n: colors_bar.append(C['blue'])
    else: colors_bar.append(C['gray'])
ax.barh(range(len(names_bar)), aucs_bar, color=colors_bar, edgecolor='black', linewidth=0.5, height=0.65)
ax.set_yticks(range(len(names_bar))); ax.set_yticklabels(names_bar, fontsize=9)
ax.set_xlabel('AUC'); ax.set_xlim(0.3, 1.0)
ax.axvline(0.5, color='gray', ls='--', lw=0.8)
for i, v in enumerate(aucs_bar):
    ax.text(v+0.01, i, f'{v:.3f}', va='center', fontsize=8)
auc_c = auc_dict.get('NetImm+MSI|EBV',0)
ax.text(0.98, -0.5, f'NetImm+MSI|EBV: {auc_c:.3f}', ha='right', fontsize=9,
        fontweight='bold', color=C['red'], transform=ax.get_xaxis_transform())
ax.set_title('C  AUC Comparison', fontsize=11, fontweight='bold', loc='left')

plt.savefig('Fig_benchmark_main.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Fig_benchmark_main.png', dpi=200, bbox_inches='tight')
plt.close(); print("✓ Fig_benchmark_main.pdf/png")

# ============ FIG 2: Heatmap ============
fig, ax = plt.subplots(figsize=(14, 5))
order_idx = scores.sort_values(['Responder','NetImmune_z'], ascending=[False,False]).index
hm_cols = ['NetImmune_z','Barrier_z','GEP18','IFNG_score','IPS','PDL1_expr']
hm_data = scores.loc[order_idx, hm_cols].copy()
for c in hm_data.columns:
    s = hm_data[c]; hm_data[c] = (s - s.mean()) / (s.std() + 1e-9)
im = ax.imshow(hm_data.values.T, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
ax.set_yticks(range(len(hm_cols)))
ax.set_yticklabels(['NetImmune','Barrier','GEP-18','IFNG','IPS','PD-L1'], fontsize=10)
ax.set_xlabel('Patients (R → NR)')
for i, r in enumerate(scores.loc[order_idx, 'Response']):
    ax.plot(i, -0.8, 's', color=resp_colors.get(r,'gray'), markersize=4, clip_on=False)
plt.colorbar(im, ax=ax, label='z-score', shrink=0.8)
ax.set_title('Score Heatmap (patients ordered by response)', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('Fig_heatmap.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Fig_heatmap.png', dpi=200, bbox_inches='tight')
plt.close(); print("✓ Fig_heatmap.pdf/png")

# ============ FIG 3: MSS subgroup ============
df_mss = scores[(~scores['MSI_H']) & (~scores['EBV_pos'])].copy()
y_mss = df_mss['Responder'].values
fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
ax = axes[0]
for i, (label, col) in enumerate([('NetImmune','NetImmune_z'),('GEP-18','GEP18'),
                                   ('IFNG','IFNG_score'),('IPS','IPS')]):
    nr_v = df_mss.loc[df_mss['Responder']==0, col].values
    r_v = df_mss.loc[df_mss['Responder']==1, col].values
    bp = ax.boxplot([nr_v], positions=[i*2.5], widths=0.6, patch_artist=True,
                    medianprops=dict(color='black'))
    bp['boxes'][0].set_facecolor('#D5D8DC')
    ax.scatter(np.random.normal(i*2.5, 0.05, len(nr_v)), nr_v, c='black', s=10, alpha=0.4, zorder=5)
    if len(r_v) > 0:
        ax.scatter([i*2.5+0.8]*len(r_v), r_v, c=C['red'], s=50, marker='*', zorder=5,
                  label='Resp' if i==0 else '')
ax.set_xticks([i*2.5 for i in range(4)])
ax.set_xticklabels(['NetImmune','GEP-18','IFNG','IPS'], fontsize=9)
ax.set_title(f'MSS/EBV- (n={len(df_mss)}, R={y_mss.sum()})', fontsize=11, fontweight='bold')
if y_mss.sum()>0: ax.legend(fontsize=8)

ax = axes[1]
if y_mss.sum() >= 2:
    for label, col, color in [('NetImmune','NetImmune_z',C['red']),('GEP-18','GEP18',C['purple']),
                               ('IFNG','IFNG_score',C['orange']),('IPS','IPS',C['green'])]:
        auc = roc_auc_score(y_mss, df_mss[col].values)
        fpr, tpr, _ = roc_curve(y_mss, df_mss[col].values)
        ax.plot(fpr, tpr, color=color, lw=2, label=f'{label} ({auc:.3f})')
    ax.plot([0,1],[0,1],'k--',lw=0.8); ax.legend(fontsize=8)
    ax.set_xlabel('FPR'); ax.set_ylabel('TPR')
    ax.set_title('ROC: MSS/EBV-', fontsize=11, fontweight='bold')
else:
    ax.text(0.5, 0.5, f'{y_mss.sum()} MSS responder(s)\nInsufficient for ROC',
            ha='center', va='center', fontsize=12, transform=ax.transAxes)
plt.tight_layout()
plt.savefig('Fig_MSS_subgroup.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Fig_MSS_subgroup.png', dpi=200, bbox_inches='tight')
plt.close(); print("✓ Fig_MSS_subgroup.pdf/png")

# =================================================================
# SUMMARY
# =================================================================
print(f"\n{'='*58}")
print(f"  FINAL RESULTS (barriernet official geneset)")
print(f"{'='*58}")
print(f"  N={len(scores)}, R={y.sum()}, NR={(1-y).sum()}")
print(f"  MSI-H={MSI_H.sum()}, EBV+={EBV_pos.sum()}\n")
for name in ['NetImm+MSI|EBV','MSI-H|EBV+','IFNG','GEP-18','PD-L1',
             'NetImmune','MSI-H','TIDE(neg)','IPS','Barrier(neg)']:
    if name in auc_dict:
        tag = ' ← OURS' if 'Net' in name or 'Barrier' in name else ''
        print(f"  {name:<20} AUC = {auc_dict[name]:.3f}{tag}")
print(f"\n  Output: benchmark_scores.tsv + 3 PDFs")
print(f"✅ Done!")
