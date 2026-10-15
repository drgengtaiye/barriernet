"""Per-line Barrier score and ln(IC50) for dasatinib and sorafenib in GDSC2 gastric lines,
reproducing GC_TME_Reprogrammer/validation/drug_vali/gdsc_validation_fix.py (same files,
same scoring), so Fig. 8g,h plot exactly the values behind C1_gdsc_fdr.csv."""
import pandas as pd, numpy as np
from scipy import stats
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation/drug_vali"
MOD = {'CAF_ECM': ['COL1A1','COL1A2','LUM','DCN','POSTN','TAGLN','ACTA2','THY1','FAP','FN1','COL3A1'],
       'TAM_SUP': ['SPP1','APOE','TREM2','C1QA','C1QB','C1QC','LGALS3','MARCO','IL10','SOCS3','STAT3']}
def zs(x):
    mu, sd = np.nanmean(x.values), np.nanstd(x.values, ddof=0)
    return x * 0.0 if (not np.isfinite(sd) or sd == 0) else (x - mu) / (sd + 1e-9)
model = pd.read_csv(f"{V}/Model.csv")
n2m = {}
for _, r in model.iterrows():
    n2m[str(r['StrippedCellLineName']).upper()] = r['ModelID']
    n2m[str(r['CellLineName']).upper().replace('-', '').replace(' ', '').replace(':', '')] = r['ModelID']
gc = set(model.loc[model['OncotreePrimaryDisease'].astype(str).str.contains('Gastric|Stomach', case=False, na=False), 'ModelID'])
ex = pd.read_csv(f"{V}/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
if 'IsDefaultEntryForModel' in ex.columns:
    ex = ex[ex['IsDefaultEntryForModel'].astype(str).str.upper().isin(['TRUE', 'YES', '1'])]
ex = ex.set_index('ModelID').select_dtypes(include=[np.number])
ex.columns = [str(c).split(' (')[0].strip().upper() for c in ex.columns]
resp = pd.read_excel(f"{V}/GDSC2_fitted_dose_response_27Oct23.xlsx")
clean = lambda s: str(s).upper().replace('-', '').replace('_', '').replace(' ', '').replace(':', '')
resp['model_id'] = resp['CELL_LINE_NAME'].apply(lambda x: n2m.get(clean(x), n2m.get(str(x).upper())))
resp = resp[resp['model_id'].isin(gc)]
common = sorted(set(ex.index) & set(resp['model_id'].unique()))
Xz = ex.loc[common].apply(zs, axis=0)
mods = pd.DataFrame({m: Xz[[g for g in G if g in Xz.columns]].mean(1) for m, G in MOD.items()})
barrier = zs(mods.mean(1))
rows = []
for drug in ["Dasatinib", "Sorafenib"]:
    ic = resp[resp['DRUG_NAME'] == drug].groupby('model_id')['LN_IC50'].median()
    c = [m for m in common if m in ic.index]
    r, p = stats.spearmanr(barrier[c], ic[c])
    print(drug, len(c), round(r, 3), round(p, 4))
    rows += [dict(drug=drug, model=m, line=model.set_index('ModelID').loc[m, 'StrippedCellLineName'],
                  barrier_z=barrier[m], ln_ic50=ic[m]) for m in c]
pd.DataFrame(rows).to_csv("Fig8_gdsc_points.csv", index=False)
