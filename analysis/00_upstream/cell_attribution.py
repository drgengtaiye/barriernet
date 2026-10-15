#!/usr/bin/env python3
"""
GSE183904 Cell-Level Attribution Analysis
==========================================
Which cell types drive Barrier vs NetImmune?

python cell_attribution.py
pip install scanpy anndata matplotlib pandas numpy scipy
==========================================
"""
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import warnings; warnings.filterwarnings('ignore')

matplotlib.rcParams.update({
    'pdf.fonttype':42,'font.family':'Arial','font.size':7,
    'axes.titlesize':8,'axes.labelsize':7,'axes.linewidth':0.5,
    'xtick.labelsize':6,'ytick.labelsize':6.5,'legend.fontsize':6,
    'figure.dpi':300,
})
mm = 1/25.4

# ======================= CONFIG =======================
H5AD_PATH = '/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/out_scanpy_strict_GSE183904/GSE183904_TumorOnly_scanpy_strict_final_scored.h5ad'
# ======================================================

# barriernet geneset_rna v1
MODULES = {
    'CAF_ECM': ['COL1A1','COL1A2','LUM','DCN','POSTN','TAGLN','ACTA2','THY1','FAP','FN1','COL3A1'],
    'TAM_SUP': ['SPP1','APOE','TREM2','C1QA','C1QB','C1QC','LGALS3','MARCO','IL10','SOCS3','STAT3'],
    'APM':     ['B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2','PSMB8','PSMB9'],
    'IFN':     ['STAT1','IRF1','CXCL9','CXCL10','IFNG'],
    'TLS_B':   ['CXCL13','LTB','MS4A1','CD79A','CD74','MZB1'],
    'T_CELL':  ['CD3D','CD3E','TRAC','TRBC1','CD8A','GZMB','NKG7'],
}

# Module → axis mapping
BARRIER_MODULES = ['CAF_ECM','TAM_SUP']
IMMUNE_MODULES  = ['APM','IFN','TLS_B','T_CELL']

# Colors for cell types (will be auto-assigned)
CMAP_CELLTYPE = {}

print("="*60)
print("  GSE183904 Cell-Level Attribution")
print("="*60)

# ============================================================
# 1. Load data
# ============================================================
print("\n--- Loading h5ad ---")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  Shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"  obs columns: {list(adata.obs.columns)}")

ct_col = "anno_level1"  # HARDCODED
ct_col = 'anno_level1'
ct_candidates = ['cell_type','celltype','CellType','cell_type_fine','cell_type_coarse',
                  'annotation','Annotation','leiden','louvain','cluster','clusters',
                  'cell_annotation','celltype_major','celltype_minor','cell_label',
                  'predicted_celltype','assigned_cell_type']
for c in ct_candidates:
    if c in adata.obs.columns:
        ct_col = c
        break
# If still not found, pick first categorical/string column that looks like annotation
if ct_col is None:
    for c in adata.obs.columns:
        if adata.obs[c].dtype == 'category' or adata.obs[c].dtype == 'object':
            nuniq = adata.obs[c].nunique()
            if 3 <= nuniq <= 50:
                ct_col = c
                break

if ct_col is None:
    print("\n  ⚠ Cannot auto-detect cell type column!")
    print("  Available columns with unique counts:")
    for c in adata.obs.columns:
        print(f"    {c}: {adata.obs[c].nunique()} unique")
    import sys; sys.exit(1)

print(f"\n  Using cell type column: '{ct_col}'")
print(f"  Cell types ({adata.obs[ct_col].nunique()}):")
for ct, n in adata.obs[ct_col].value_counts().items():
    print(f"    {ct}: {n} cells ({n/len(adata)*100:.1f}%)")

# ============================================================
# 2. Check expression layer
# ============================================================
print("\n--- Expression data ---")

# Use .raw if available (usually has all genes before HVG filtering)
if adata.raw is not None:
    print("  Using adata.raw (pre-HVG)")
    expr_adata = adata.raw.to_adata()
else:
    expr_adata = adata

# Get gene names
gene_names = [str(g).upper() for g in expr_adata.var_names]
expr_adata.var_names = gene_names

# Check key genes
all_module_genes = set()
for genes in MODULES.values():
    all_module_genes.update([g.upper() for g in genes])
found = all_module_genes & set(gene_names)
print(f"  Module genes found: {len(found)}/{len(all_module_genes)}")

missing = all_module_genes - set(gene_names)
if missing:
    print(f"  Missing: {missing}")

# ============================================================
# 3. Compute module scores per cell
# ============================================================
print("\n--- Computing module scores ---")

# Get expression matrix (dense for small gene sets)
from scipy.sparse import issparse

def get_expr(adata_obj, genes):
    """Get expression values for genes, handling sparse matrices."""
    avail = [g for g in genes if g in adata_obj.var_names]
    if len(avail) == 0:
        return None, avail
    idx = [list(adata_obj.var_names).index(g) for g in avail]
    X = adata_obj.X[:, idx]
    if issparse(X):
        X = X.toarray()
    return X, avail

for m, genes in MODULES.items():
    genes_upper = [g.upper() for g in genes]
    X, avail = get_expr(expr_adata, genes_upper)
    if X is None or len(avail) < 3:
        print(f"  {m}: {len(avail) if avail else 0}/{len(genes)} genes ⚠ SKIP")
        adata.obs[m] = np.nan
        continue
    # Score = mean expression (on log-normalized scale)
    adata.obs[m] = np.mean(X, axis=1)
    print(f"  {m}: {len(avail)}/{len(genes)} genes ✓")

# Barrier / Immune / NetImmune per cell
barrier_cols = [m for m in BARRIER_MODULES if m in adata.obs.columns and adata.obs[m].notna().any()]
immune_cols  = [m for m in IMMUNE_MODULES if m in adata.obs.columns and adata.obs[m].notna().any()]

adata.obs['Barrier']   = adata.obs[barrier_cols].mean(axis=1)
adata.obs['Immune']    = adata.obs[immune_cols].mean(axis=1)
adata.obs['NetImmune'] = adata.obs['Immune'] - adata.obs['Barrier']

# ============================================================
# 4. Cell-type × module attribution table
# ============================================================
print("\n--- Attribution analysis ---")

cell_types = adata.obs[ct_col].unique()
score_cols = list(MODULES.keys()) + ['Barrier','Immune','NetImmune']

# Mean score per cell type
attr_df = pd.DataFrame(index=cell_types)
for col in score_cols:
    if col in adata.obs.columns:
        attr_df[col] = adata.obs.groupby(ct_col)[col].mean()

# Z-score across cell types (for heatmap)
attr_z = attr_df.copy()
for col in attr_z.columns:
    mu, sd = attr_z[col].mean(), attr_z[col].std()
    attr_z[col] = (attr_z[col] - mu) / (sd + 1e-9) if sd > 0 else 0

# Sort cell types by NetImmune (low→high = Barrier-dominant → Immune-dominant)
attr_z = attr_z.sort_values('NetImmune')
attr_df = attr_df.loc[attr_z.index]

# Print
print(f"\n  {'Cell Type':<25} {'CAF_ECM':>8} {'TAM_SUP':>8} {'APM':>8} {'IFN':>8} {'TLS_B':>8} {'T_CELL':>8} │ {'Barrier':>8} {'NetImm':>8}")
print(f"  {'-'*105}")
for ct in attr_z.index:
    row = attr_df.loc[ct]
    print(f"  {str(ct):<25} {row.get('CAF_ECM',0):>8.3f} {row.get('TAM_SUP',0):>8.3f} {row.get('APM',0):>8.3f} {row.get('IFN',0):>8.3f} {row.get('TLS_B',0):>8.3f} {row.get('T_CELL',0):>8.3f} │ {row.get('Barrier',0):>8.3f} {row.get('NetImmune',0):>8.3f}")

attr_df.to_csv('cell_attribution_scores.csv')
attr_z.to_csv('cell_attribution_zscore.csv')
print("\n  ✓ Saved: cell_attribution_scores.csv, cell_attribution_zscore.csv")

# ============================================================
# 5. Fraction-weighted attribution
# ============================================================
print("\n--- Fraction-weighted contribution ---")
ct_counts = adata.obs[ct_col].value_counts()
ct_frac = ct_counts / ct_counts.sum()

# Contribution = fraction × mean score
contrib = pd.DataFrame(index=attr_df.index)
for col in MODULES.keys():
    if col in attr_df.columns:
        contrib[col] = ct_frac[attr_df.index] * attr_df[col]

# Normalize per module
contrib_norm = contrib.copy()
for col in contrib_norm.columns:
    total = contrib_norm[col].abs().sum()
    if total > 0:
        contrib_norm[col] = contrib_norm[col] / total * 100

print(f"\n  Contribution (%) per module:")
print(f"  {'Cell Type':<25} {'CAF_ECM':>8} {'TAM_SUP':>8} {'APM':>8} {'IFN':>8} {'TLS_B':>8} {'T_CELL':>8}")
print(f"  {'-'*80}")
for ct in contrib_norm.index:
    row = contrib_norm.loc[ct]
    print(f"  {str(ct):<25} {row.get('CAF_ECM',0):>7.1f}% {row.get('TAM_SUP',0):>7.1f}% {row.get('APM',0):>7.1f}% {row.get('IFN',0):>7.1f}% {row.get('TLS_B',0):>7.1f}% {row.get('T_CELL',0):>7.1f}%")

contrib_norm.to_csv('cell_attribution_contribution_pct.csv')

# ============================================================
# 6. FIGURES
# ============================================================
print("\n--- Generating figures ---")

# Assign colors to cell types
n_ct = len(attr_z.index)
if n_ct <= 10:
    cpal = plt.cm.tab10(np.linspace(0, 1, 10))[:n_ct]
elif n_ct <= 20:
    cpal = plt.cm.tab20(np.linspace(0, 1, 20))[:n_ct]
else:
    cpal = plt.cm.gist_ncar(np.linspace(0.05, 0.95, n_ct))
ct_colors = {ct: cpal[i] for i, ct in enumerate(attr_z.index)}

# ===== Figure: 180mm × 160mm, 2×2 panels =====
fig = plt.figure(figsize=(180*mm, 160*mm))
gs = gridspec.GridSpec(2, 2, hspace=0.55, wspace=0.45)

# ----- a: Heatmap (cell type × module z-scores) -----
ax = fig.add_subplot(gs[0, :])

hm_cols = ['CAF_ECM','TAM_SUP','APM','IFN','TLS_B','T_CELL']
hm_data = attr_z[hm_cols].values

im = ax.imshow(hm_data, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
ax.set_yticks(range(n_ct))
ax.set_yticklabels(attr_z.index, fontsize=5.5)
ax.set_xticks(range(len(hm_cols)))
ax.set_xticklabels(hm_cols, fontsize=6.5, rotation=30, ha='right')

# Add value text
for i in range(hm_data.shape[0]):
    for j in range(hm_data.shape[1]):
        val = hm_data[i, j]
        color = 'white' if abs(val) > 1.2 else 'black'
        ax.text(j, i, f'{val:.1f}', ha='center', va='center', fontsize=4.5, color=color)

# Barrier / Immune brackets
ax.annotate('', xy=(-0.5, -0.8), xytext=(1.5, -0.8),
            arrowprops=dict(arrowstyle='-', color='#C0392B', lw=1.5), annotation_clip=False)
ax.text(0.5, -1.2, 'Barrier', ha='center', fontsize=6, fontweight='bold', color='#C0392B',
        clip_on=False)
ax.annotate('', xy=(2, -0.8), xytext=(5.5, -0.8),
            arrowprops=dict(arrowstyle='-', color='#2980B9', lw=1.5), annotation_clip=False)
ax.text(3.75, -1.2, 'Immune', ha='center', fontsize=6, fontweight='bold', color='#2980B9',
        clip_on=False)

cb = plt.colorbar(im, ax=ax, shrink=0.6, aspect=12, pad=0.02)
cb.ax.tick_params(labelsize=5); cb.set_label('z-score', fontsize=6)
ax.set_title('Module enrichment across cell types', fontweight='bold', pad=12)
ax.text(-0.08, 1.08, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold')
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

# ----- b: Stacked bar (contribution %) -----
ax = fig.add_subplot(gs[1, 0])

# Barrier contribution
barrier_contrib = contrib_norm[['CAF_ECM','TAM_SUP']].copy()
barrier_contrib['Barrier_total'] = barrier_contrib.sum(axis=1)
barrier_top = barrier_contrib.nlargest(8, 'Barrier_total')

# Plot stacked
bottom = np.zeros(len(barrier_top))
bars_ecm = ax.barh(range(len(barrier_top)), barrier_top['CAF_ECM'], height=0.6,
                    color='#E74C3C', alpha=0.7, label='CAF_ECM')
bottom += barrier_top['CAF_ECM'].values
bars_tam = ax.barh(range(len(barrier_top)), barrier_top['TAM_SUP'], left=bottom, height=0.6,
                    color='#E67E22', alpha=0.7, label='TAM_SUP')

ax.set_yticks(range(len(barrier_top)))
ax.set_yticklabels(barrier_top.index, fontsize=6)
ax.set_xlabel('Contribution (%)', fontsize=7)
ax.set_title('Top contributors to Barrier', fontweight='bold', pad=4, color='#C0392B')
ax.legend(fontsize=5.5, frameon=False)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.text(-0.2, 1.05, 'b', transform=ax.transAxes, fontsize=10, fontweight='bold')

# ----- c: Stacked bar (Immune contribution) -----
ax = fig.add_subplot(gs[1, 1])

immune_contrib = contrib_norm[['APM','IFN','TLS_B','T_CELL']].copy()
immune_contrib['Immune_total'] = immune_contrib.sum(axis=1)
immune_top = immune_contrib.nlargest(8, 'Immune_total')

bottom = np.zeros(len(immune_top))
imm_colors = {'APM':'#3498DB','IFN':'#1ABC9C','TLS_B':'#9B59B6','T_CELL':'#2ECC71'}
for col in ['APM','IFN','TLS_B','T_CELL']:
    ax.barh(range(len(immune_top)), immune_top[col], left=bottom, height=0.6,
            color=imm_colors[col], alpha=0.7, label=col)
    bottom += immune_top[col].values

ax.set_yticks(range(len(immune_top)))
ax.set_yticklabels(immune_top.index, fontsize=6)
ax.set_xlabel('Contribution (%)', fontsize=7)
ax.set_title('Top contributors to Immune', fontweight='bold', pad=4, color='#2980B9')
ax.legend(fontsize=5.5, frameon=False, ncol=2)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.text(-0.2, 1.05, 'c', transform=ax.transAxes, fontsize=10, fontweight='bold')

plt.savefig('Fig_cell_attribution.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Fig_cell_attribution.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Fig_cell_attribution.pdf/png")

# ===== Individual panels =====
# Panel a alone
fig, ax = plt.subplots(figsize=(160*mm, max(40*mm, n_ct*4*mm)))
im = ax.imshow(hm_data, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2, interpolation='nearest')
ax.set_yticks(range(n_ct)); ax.set_yticklabels(attr_z.index, fontsize=5.5)
ax.set_xticks(range(len(hm_cols))); ax.set_xticklabels(hm_cols, fontsize=7, rotation=30, ha='right')
for i in range(hm_data.shape[0]):
    for j in range(hm_data.shape[1]):
        val = hm_data[i,j]; color = 'white' if abs(val)>1.2 else 'black'
        ax.text(j, i, f'{val:.1f}', ha='center', va='center', fontsize=5, color=color)
cb = plt.colorbar(im, ax=ax, shrink=0.6, aspect=12, pad=0.02)
cb.ax.tick_params(labelsize=5); cb.set_label('z-score', fontsize=6)
ax.set_title('Module enrichment across cell types', fontweight='bold', pad=6)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.savefig('panel_attribution_heatmap.pdf', dpi=300, bbox_inches='tight')
plt.savefig('panel_attribution_heatmap.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ panel_attribution_heatmap.pdf/png")

# ===== Violin plots for key modules in key cell types =====
fig, axes = plt.subplots(1, 3, figsize=(170*mm, 55*mm))

violin_data = [
    ('CAF_ECM', 'Barrier component', '#E74C3C'),
    ('TAM_SUP', 'Barrier component', '#E67E22'),
    ('T_CELL',  'Immune component',  '#2ECC71'),
]

for idx, (mod, subtitle, color) in enumerate(violin_data):
    ax = axes[idx]
    if mod not in adata.obs.columns:
        continue
    
    # Top 6 cell types by this module score
    ct_means = adata.obs.groupby(ct_col)[mod].mean().sort_values(ascending=False)
    top_cts = ct_means.head(6).index.tolist()
    
    data_list = []
    labels = []
    for ct in top_cts:
        vals = adata.obs.loc[adata.obs[ct_col]==ct, mod].values
        data_list.append(vals)
        labels.append(str(ct)[:15])  # truncate long names
    
    parts = ax.violinplot(data_list, positions=range(len(data_list)),
                          showmeans=True, showmedians=False, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(color); pc.set_alpha(0.5)
    parts['cmeans'].set_color('black'); parts['cmeans'].set_linewidth(0.8)
    
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=5, rotation=40, ha='right')
    ax.set_ylabel('Score', fontsize=7)
    ax.set_title(mod, fontweight='bold', fontsize=8, pad=4)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('panel_attribution_violin.pdf', dpi=300, bbox_inches='tight')
plt.savefig('panel_attribution_violin.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ panel_attribution_violin.pdf/png")

# ============================================================
# 7. Summary
# ============================================================
print(f"\n{'='*60}")
print(f"  SUMMARY")
print(f"{'='*60}")

# Top Barrier drivers
barrier_mean = attr_df[['CAF_ECM','TAM_SUP']].mean(axis=1).sort_values(ascending=False)
print(f"\n  Top Barrier-driving cell types:")
for ct in barrier_mean.head(3).index:
    print(f"    {ct}: CAF_ECM={attr_df.loc[ct,'CAF_ECM']:.3f}, TAM_SUP={attr_df.loc[ct,'TAM_SUP']:.3f}")

# Top Immune drivers
immune_mean = attr_df[['APM','IFN','TLS_B','T_CELL']].mean(axis=1).sort_values(ascending=False)
print(f"\n  Top Immune-driving cell types:")
for ct in immune_mean.head(3).index:
    r = attr_df.loc[ct]
    print(f"    {ct}: APM={r['APM']:.3f}, IFN={r['IFN']:.3f}, TLS_B={r['TLS_B']:.3f}, T_CELL={r['T_CELL']:.3f}")

print(f"\n  Output: 3 CSVs + Fig_cell_attribution.pdf + panel_*.pdf")
print(f"✅ Done!")
