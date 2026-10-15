"""NDCG@5 of the leave-one-section-out ranking model (Fig7_AI_rank) against a leave-one-
section-out global-mean baseline (each compound's mean score in the other nine sections).
Relevance = true drug score minus the section minimum; standard log2 discount."""
import numpy as np, pandas as pd
R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
S = pd.read_csv(f"{R}/results/ai_recommendations_drug_top15_clean.csv", index_col=0)
B = pd.read_csv(f"{R}/results/ring_barrier_state_v2/BarrierIndex_by_slice_with_clinical.csv")
S.index = S.index.map(dict(zip(B["slice"], B["GC_raw"])))
P = pd.read_csv(f"{R}/figures/Fig7_AI_rank/Fig7_rank_pred_scores_LOO.csv", index_col=0)[S.columns]
def ndcg(true, pred, k=5):
    rel = true - true.min()
    order = pred.sort_values(ascending=False).index[:k]
    ideal = rel.sort_values(ascending=False).index[:k]
    disc = 1 / np.log2(np.arange(2, k + 2))
    return (rel[order].values * disc).sum() / (rel[ideal].values * disc).sum()
rows = []
for s in S.index:
    base = S.drop(index=s).mean()
    rows.append(dict(section=s, model=ndcg(S.loc[s], P.loc[s]), baseline=ndcg(S.loc[s], base),
                     top5_model=len(set(P.loc[s].nlargest(5).index) & set(S.loc[s].nlargest(5).index)) / 5,
                     top5_baseline=len(set(base.nlargest(5).index) & set(S.loc[s].nlargest(5).index)) / 5))
D = pd.DataFrame(rows); D.to_csv("Fig8_rank_vs_baseline.csv", index=False)
print(D.round(3).to_string(index=False)); print(D.mean(numeric_only=True).round(3))
