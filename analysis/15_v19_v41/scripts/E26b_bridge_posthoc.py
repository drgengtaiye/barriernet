"""E26b (post hoc, not pre-registered) - why did E26 not replicate the published Visium reading?
A  published sections (GSE251950, 10) with the published inputs: Tangram states, TLS_niche_score peaks
B  published sections with the E26 inputs: barriernet module z-scores, TLS_B peaks
C  E26 BarrierIndex against stage at donor level (stage is a donor property; 20 donors)
Rings: barriernet 0.2.1 assign_rings defaults (kNN-9 peaks at the 95th percentile, Far >= 3 hops)."""
import warnings, glob, sys; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import binomtest, spearmanr
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/repo_barriernet/barriernet")
from barriernet.spatial import assign_rings
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2"
MOD = dict(CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
           TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
           SPP1_arm=["SPP1", "TREM2", "MARCO", "APOE", "LGALS3"])
rows = []
for f in sorted(glob.glob(f"{V}/*.with_state_and_ring.h5ad")):
    a = ad.read_h5ad(f); obs = a.obs.copy(); xy = np.asarray(a.obsm["spatial"], float)
    need = [g for g in sum(MOD.values(), []) if g in a.var_names]
    X = pd.DataFrame(a[:, need].X.toarray() if hasattr(a.X, "toarray") else a[:, need].X, index=a.obs_names, columns=need)
    for k, gs in MOD.items():
        g = [x for x in gs if x in X.columns and X[x].std() > 0]
        obs[f"m_{k}"] = ((X[g] - X[g].mean()) / X[g].std()).mean(1) if len(g) >= 3 else np.nan
    rA = assign_rings(obs, xy, "TLS_niche_score").ring.values
    rB = assign_rings(obs, xy, "m_TLS_B").ring.values
    d = lambda v, r: np.nanmean(v[r == "Ring1"]) - np.nanmean(v[r == "Far"])
    rows.append(dict(section=f.split("/")[-1].split(".")[0], spots=a.n_obs,
                     A_state_SPP1_TAM=d(obs.state_SPP1_TAM.values, rA), A_state_CAF_ECM=d(obs.state_CAF_ECM.values, rA),
                     B_module_SPP1_arm=d(obs.m_SPP1_arm.values, rB), B_module_CAF_ECM=d(obs.m_CAF_ECM.values, rB),
                     mixed_module_SPP1_arm_on_published_peaks=d(obs.m_SPP1_arm.values, rA),
                     r_TLSniche_TLSB=np.corrcoef(obs.TLS_niche_score, obs.m_TLS_B)[0, 1],
                     r_stateSPP1_moduleSPP1=np.corrcoef(obs.state_SPP1_TAM, obs.m_SPP1_arm)[0, 1]))
B = pd.DataFrame(rows); B.to_csv(f"{O}/E26b_bridge_sections.csv", index=False)
print(B.round(3).to_string(index=False)); print()
for c in [c for c in B.columns if c.startswith(("A_", "B_", "mixed"))]:
    neg = int((B[c] < 0).sum()); n = int(B[c].notna().sum())
    print(f"{c:42s} Ring1 < Far in {neg}/{n}; sign test (two-sided) P = {binomtest(neg, n, 0.5).pvalue:.4f}")
print(f"median spot correlation: TLS_niche vs TLS_B module {B.r_TLSniche_TLSB.median():.2f}; "
      f"SPP1_TAM state vs SPP1_arm module {B.r_stateSPP1_moduleSPP1.median():.2f}")
# C: donor-level BarrierIndex vs stage in E26
T = pd.read_csv(f"{O}/E26_sections.csv")
Dn = T.groupby("donor").agg(BarrierIndex=("BarrierIndex", "mean"), stage=("stage", "first"), n=("gsm", "size"))
st = Dn.stage.map({"Stage I": 1, "Stage II": 2, "Stage III": 3})
rho, p = spearmanr(Dn.BarrierIndex, st)
print(f"\nC  E26 BarrierIndex vs stage, donor level: rho {rho:.2f}, P = {p:.3f} (20 donors; section level was rho -0.31, P = 0.010)")
print("   stage counts (donors):", Dn.stage.value_counts().to_dict())
