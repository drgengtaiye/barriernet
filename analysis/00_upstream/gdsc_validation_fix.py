#!/usr/bin/env python3
"""
GDSC Validation (fixed for DepMap 25Q3 format)
"""
import pandas as pd, numpy as np
from scipy import stats
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings; warnings.filterwarnings('ignore')

matplotlib.rcParams.update({
    'pdf.fonttype':42,'font.family':'Arial','font.size':7,
    'axes.titlesize':8,'axes.labelsize':7,'axes.linewidth':0.5,
    'figure.dpi':300,
})
mm = 1/25.4

MODULES = {
    'CAF_ECM':['COL1A1','COL1A2','LUM','DCN','POSTN','TAGLN','ACTA2','THY1','FAP','FN1','COL3A1'],
    'TAM_SUP':['SPP1','APOE','TREM2','C1QA','C1QB','C1QC','LGALS3','MARCO','IL10','SOCS3','STAT3'],
    'APM':['B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2','PSMB8','PSMB9'],
    'IFN':['STAT1','IRF1','CXCL9','CXCL10','IFNG'],
    'TLS_B':['CXCL13','LTB','MS4A1','CD79A','CD74','MZB1'],
    'T_CELL':['CD3D','CD3E','TRAC','TRBC1','CD8A','GZMB','NKG7'],
}

TOP_DRUGS = {
    'Dasatinib':['Dasatinib'], 'SB431542':['SB431542'],
    'Trametinib':['Trametinib'], 'Cediranib':['Cediranib','AZD2171'],
    'XAV-939':['XAV939','XAV-939'], 'Nilotinib':['Nilotinib'],
    'Bosutinib':['Bosutinib'], 'Vorinostat':['Vorinostat'],
    'Sorafenib':['Sorafenib'], 'Erlotinib':['Erlotinib'],
    'Lapatinib':['Lapatinib'], 'PD-0325901':['PD0325901','PD-0325901'],
    'Imatinib':['Imatinib'], 'Enzastaurin':['Enzastaurin'],
    'Thalidomide':['Thalidomide'], 'Maraviroc':['Maraviroc'],
    'Curcumin':['Curcumin'],
}

def zs(x):
    x=pd.to_numeric(x,errors='coerce')
    mu,sd=np.nanmean(x.values),np.nanstd(x.values,ddof=0)
    return x*0.0 if (not np.isfinite(sd) or sd==0) else (x-mu)/(sd+1e-9)

print("="*65)
print("  GDSC Drug Sensitivity Validation")
print("="*65)

# ============================================================
# 1. Load Model.csv → build ModelID ↔ CellLineName mapping
# ============================================================
print("\n--- Loading Model.csv ---")
model = pd.read_csv('Model.csv')
# Build mapping: ModelID → StrippedCellLineName (uppercase, no special chars)
model_map = dict(zip(model['ModelID'], model['StrippedCellLineName']))
# Also: StrippedCellLineName → ModelID (reverse)
name_to_model = {}
for _, row in model.iterrows():
    name_to_model[str(row['StrippedCellLineName']).upper()] = row['ModelID']
    name_to_model[str(row['CellLineName']).upper().replace('-','').replace(' ','').replace(':','')] = row['ModelID']

# Get gastric cell lines
gc_models = model[model['OncotreePrimaryDisease'].astype(str).str.contains(
    'Gastric|Stomach', case=False, na=False)]
print(f"  Total models: {len(model)}")
print(f"  Gastric models: {len(gc_models)}")
if len(gc_models) > 0:
    print(f"  GC cell lines: {list(gc_models['CellLineName'].head(10))}")
gc_model_ids = set(gc_models['ModelID'].values)

# ============================================================
# 2. Load expression (DepMap format)
# ============================================================
print("\n--- Loading expression ---")
expr_raw = pd.read_csv('OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv')
print(f"  Raw shape: {expr_raw.shape}")

# Use ModelID as index, drop metadata columns
meta_cols = ['Unnamed: 0','SequencingID','ModelID','IsDefaultEntryForModel',
             'ModelConditionID','IsDefaultEntryForMC','ProfileID']
id_col = 'ModelID'

# Keep only default entries per model
if 'IsDefaultEntryForModel' in expr_raw.columns:
    expr_raw = expr_raw[expr_raw['IsDefaultEntryForModel'].astype(str).str.upper().isin(['TRUE','YES','1'])]
    print(f"  After default filter: {expr_raw.shape[0]} models")

expr = expr_raw.set_index(id_col)
drop_cols = [c for c in meta_cols if c in expr.columns and c != id_col]
expr = expr.drop(columns=drop_cols, errors='ignore')
expr = expr.select_dtypes(include=[np.number])

# Clean gene names: "GENE (ENTREZ)" → "GENE"
expr.columns = [str(c).split(' (')[0].strip().upper() for c in expr.columns]
print(f"  Expression: {expr.shape[0]} models × {expr.shape[1]} genes")

key_genes = ['COL1A1','SPP1','B2M','STAT1','CXCL13','CD8A']
found = [g for g in key_genes if g in expr.columns]
print(f"  Key genes: {len(found)}/6 ({found})")

# ============================================================
# 3. Load GDSC2 drug response
# ============================================================
print("\n--- Loading GDSC2 ---")
resp = pd.read_excel('GDSC2_fitted_dose_response_27Oct23.xlsx')
print(f"  Shape: {resp.shape}")

# Map GDSC cell line names to DepMap ModelIDs
def clean(s):
    return str(s).upper().replace('-','').replace('_','').replace(' ','').replace(':','')

resp['model_id'] = resp['CELL_LINE_NAME'].apply(
    lambda x: name_to_model.get(clean(x), name_to_model.get(str(x).upper(), None)))

matched_resp = resp[resp['model_id'].notna()].copy()
print(f"  Matched to DepMap: {matched_resp['model_id'].nunique()} cell lines "
      f"({len(matched_resp)} entries)")

# Filter gastric
gc_resp = matched_resp[matched_resp['model_id'].isin(gc_model_ids)]
print(f"  Gastric entries: {len(gc_resp)} ({gc_resp['model_id'].nunique()} cell lines)")

# If too few gastric, use pan-cancer
if gc_resp['model_id'].nunique() < 10:
    print(f"  Too few GC lines → using pan-cancer ({matched_resp['model_id'].nunique()} lines)")
    use_resp = matched_resp
    analysis_scope = 'pan-cancer'
else:
    use_resp = gc_resp
    analysis_scope = 'gastric'

# Match drugs
available_drugs = set(use_resp['DRUG_NAME'].unique())
matched_drugs = {}
for our_name, aliases in TOP_DRUGS.items():
    for alias in aliases:
        hits = [d for d in available_drugs if alias.lower() == d.lower()]
        if hits:
            matched_drugs[our_name] = hits[0]
            break
        hits = [d for d in available_drugs if alias.lower() in d.lower()]
        if hits:
            matched_drugs[our_name] = hits[0]
            break

print(f"\n  Matched drugs: {len(matched_drugs)}")
for k, v in matched_drugs.items():
    print(f"    {k} → {v}")

# ============================================================
# 4. Compute barriernet scores
# ============================================================
print("\n--- barriernet scores ---")
common_models = sorted(set(expr.index) & set(use_resp['model_id'].unique()))
print(f"  Models with both expr + IC50: {len(common_models)}")

expr_sub = expr.loc[[m for m in common_models if m in expr.index]]
Xz = expr_sub.apply(zs, axis=0)

def zmod(genes):
    a = [g.upper() for g in genes if g.upper() in Xz.columns]
    return Xz[a].mean(axis=1) if len(a)>=3 else pd.Series(np.nan, index=Xz.index)

mod_df = pd.DataFrame({m: zmod(g) for m,g in MODULES.items()})
mod_df['Barrier'] = mod_df[['CAF_ECM','TAM_SUP']].mean(axis=1)
mod_df['Immune'] = mod_df[['APM','IFN','TLS_B','T_CELL']].mean(axis=1)
mod_df['NetImmune'] = mod_df['Immune'] - mod_df['Barrier']
mod_df['Barrier_z'] = zs(mod_df['Barrier'])
mod_df['NetImmune_z'] = zs(mod_df['NetImmune'])

for m in MODULES:
    avail = [g.upper() for g in MODULES[m] if g.upper() in Xz.columns]
    print(f"  {m}: {len(avail)}/{len(MODULES[m])}")
print(f"  Barrier range: [{mod_df['Barrier'].min():.2f}, {mod_df['Barrier'].max():.2f}]")

# ============================================================
# 5. Correlate Barrier vs IC50
# ============================================================
print(f"\n{'='*60}")
print(f"  RESULTS: Barrier vs Drug Sensitivity ({analysis_scope})")
print(f"{'='*60}")
print(f"  {'Drug':<20} {'n':>4} {'Barrier_r':>10} {'p':>10} {'NetImm_r':>10} {'p':>10}")
print(f"  {'-'*66}")

results = []
for our_name, gdsc_name in matched_drugs.items():
    drug_data = use_resp[use_resp['DRUG_NAME'] == gdsc_name][['model_id','LN_IC50']].copy()
    drug_data = drug_data.groupby('model_id')['LN_IC50'].median()  # median if duplicates
    
    common = mod_df.index.intersection(drug_data.index)
    if len(common) < 8:
        continue
    
    ic50 = drug_data[common].astype(float)
    barrier = mod_df.loc[common, 'Barrier_z']
    netimm = mod_df.loc[common, 'NetImmune_z']
    
    valid = ic50.notna() & barrier.notna()
    if valid.sum() < 8:
        continue
    
    r_bar, p_bar = stats.spearmanr(barrier[valid], ic50[valid])
    r_ni, p_ni = stats.spearmanr(netimm[valid], ic50[valid])
    
    results.append({
        'Drug': our_name, 'GDSC_name': gdsc_name, 'n': valid.sum(),
        'Barrier_r': r_bar, 'Barrier_p': p_bar,
        'NetImmune_r': r_ni, 'NetImmune_p': p_ni,
    })
    
    sb = '**' if p_bar<0.01 else '*' if p_bar<0.05 else ''
    sn = '**' if p_ni<0.01 else '*' if p_ni<0.05 else ''
    hit = '← ✓' if r_bar < 0 and p_bar < 0.05 else ''
    print(f"  {our_name:<20} {valid.sum():>4} {r_bar:>+10.3f}{sb:<3} {p_bar:>9.4f} {r_ni:>+10.3f}{sn:<3} {p_ni:>9.4f}  {hit}")

res_df = pd.DataFrame(results)
res_df.to_csv('gdsc_validation_results.csv', index=False)

# Summary
n_neg = (res_df['Barrier_r'] < 0).sum()
n_sig = ((res_df['Barrier_r'] < 0) & (res_df['Barrier_p'] < 0.05)).sum()
print(f"\n  Negative Barrier-IC50 correlation: {n_neg}/{len(res_df)} ({n_neg/len(res_df)*100:.0f}%)")
print(f"  Significant (p<0.05): {n_sig}/{len(res_df)}")

# ============================================================
# 6. Figure
# ============================================================
print("\n--- Figure ---")
C_BAR = '#C0392B'

res_sorted = res_df.sort_values('Barrier_r')
n = len(res_sorted)

fig, axes = plt.subplots(1, 2, figsize=(170*mm, max(60, n*8+20)*mm),
                          gridspec_kw={'width_ratios':[1.2, 1], 'wspace':0.4})

# Panel a: forest plot
ax = axes[0]
for i, (_, row) in enumerate(res_sorted.iterrows()):
    r = row['Barrier_r']
    col = C_BAR if r < 0 else '#95A5A6'
    alpha = 0.9 if row['Barrier_p'] < 0.05 else 0.4
    sig = '**' if row['Barrier_p']<0.01 else '*' if row['Barrier_p']<0.05 else ''
    
    ax.barh(i, r, color=col, alpha=alpha, height=0.65, edgecolor='white', linewidth=0.3)
    lx = r - 0.02 if r < 0 else r + 0.02
    ha = 'right' if r < 0 else 'left'
    ax.text(lx, i, f'{r:.3f}{sig}', va='center', ha=ha, fontsize=5.5,
            fontweight='bold' if sig else 'normal')

ax.axvline(0, color='#999', lw=0.5, ls='--')
ax.set_yticks(range(n))
ax.set_yticklabels([f"{r['Drug']} (n={r['n']})" for _,r in res_sorted.iterrows()], fontsize=6)
ax.set_xlabel('Spearman r (Barrier vs ln(IC50))', fontsize=7)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.set_title(f'Barrier score vs drug sensitivity\n({analysis_scope}, negative = more sensitive)',
             fontweight='bold', pad=6)
ax.text(-0.3, 1.03, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold')

# Panel b: best scatter
ax = axes[1]
neg_sig = res_df[res_df['Barrier_r'] < 0].sort_values('Barrier_p')

if len(neg_sig) > 0:
    best = neg_sig.iloc[0]
    drug_data = use_resp[use_resp['DRUG_NAME'] == best['GDSC_name']][['model_id','LN_IC50']]
    drug_data = drug_data.groupby('model_id')['LN_IC50'].median()
    common = mod_df.index.intersection(drug_data.index)
    ic50 = drug_data[common].astype(float)
    barrier = mod_df.loc[common, 'Barrier_z']
    valid = ic50.notna() & barrier.notna()
    
    # Color by gastric vs other
    colors = []
    for m in common[valid]:
        colors.append(C_BAR if m in gc_model_ids else '#95A5A6')
    
    ax.scatter(barrier[valid], ic50[valid], c=colors, s=12, alpha=0.6,
               edgecolors='white', linewidths=0.2)
    
    slope, intercept, r, p, se = stats.linregress(barrier[valid], ic50[valid])
    xl = np.linspace(barrier[valid].min(), barrier[valid].max(), 100)
    ax.plot(xl, slope*xl+intercept, color=C_BAR, lw=1.2, ls='--')
    
    ax.set_xlabel('Barrier score (z)', fontsize=7)
    ax.set_ylabel(f'ln(IC50)', fontsize=7)
    ax.set_title(f'{best["Drug"]}\nr={best["Barrier_r"]:.3f}, p={best["Barrier_p"]:.4f}',
                 fontweight='bold', pad=6)
    
    from matplotlib.lines import Line2D
    ax.legend([Line2D([0],[0],marker='o',color='w',mfc=C_BAR,ms=4),
               Line2D([0],[0],marker='o',color='w',mfc='#95A5A6',ms=4)],
              ['Gastric','Other'], fontsize=5.5, frameon=False, loc='upper right')
else:
    ax.text(0.5, 0.5, 'No negative correlations', ha='center', va='center',
            transform=ax.transAxes, fontsize=10)

ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.text(-0.18, 1.03, 'b', transform=ax.transAxes, fontsize=10, fontweight='bold')

plt.savefig('Fig_GDSC_validation.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Fig_GDSC_validation.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Fig_GDSC_validation.pdf/png")

# Top 4 scatter
neg_top = res_df[res_df['Barrier_r'] < 0].sort_values('Barrier_p').head(4)
if len(neg_top) >= 2:
    fig, axes2 = plt.subplots(1, len(neg_top), figsize=(len(neg_top)*45*mm, 45*mm))
    if len(neg_top)==1: axes2=[axes2]
    for idx, (_, row) in enumerate(neg_top.iterrows()):
        ax = axes2[idx]
        dd = use_resp[use_resp['DRUG_NAME']==row['GDSC_name']][['model_id','LN_IC50']]
        dd = dd.groupby('model_id')['LN_IC50'].median()
        cm = mod_df.index.intersection(dd.index)
        ic = dd[cm].astype(float); br = mod_df.loc[cm,'Barrier_z']
        v = ic.notna() & br.notna()
        ax.scatter(br[v], ic[v], c=C_BAR, s=10, alpha=0.4, edgecolors='white', linewidths=0.2)
        sl,it,_,_,_ = stats.linregress(br[v],ic[v])
        xl2 = np.linspace(br[v].min(),br[v].max(),50)
        ax.plot(xl2, sl*xl2+it, color=C_BAR, lw=1, ls='--')
        sig = '**' if row['Barrier_p']<0.01 else '*' if row['Barrier_p']<0.05 else ''
        ax.set_title(f"{row['Drug']}\nr={row['Barrier_r']:.3f}{sig}", fontsize=7, fontweight='bold')
        ax.set_xlabel('Barrier (z)', fontsize=6)
        if idx==0: ax.set_ylabel('ln(IC50)', fontsize=6)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig('Fig_GDSC_scatter_top4.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('Fig_GDSC_scatter_top4.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Fig_GDSC_scatter_top4.pdf/png")

print(f"\n✅ Done! → gdsc_validation_results.csv + figures")
