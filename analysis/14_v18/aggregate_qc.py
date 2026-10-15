"""PRJEB60680 baseline: aggregate salmon transcript quantifications to genes and run QC (no outcome is used here).

Inputs : quant/salmon/<E##>/quant.sf (primary: fastp-trimmed reads; salmon 2.0.0, Ensembl 110 GRCh38 cDNA index, k = 31),
         quant/salmon_untrimmed/<E##>/quant.sf (sensitivity), quant/fastp/<run>.json, quant/E36_mate_check/,
         quant/tx2gene_Ensembl110_cDNA.tsv.gz (from the cDNA FASTA headers that built the index),
         fastq/baseline_files.tsv (ENA filereport), fastq/baseline_instrument.tsv (instrument from submitted read headers),
         ../Wright2026_CellRepMed_TableS1_clinical_47pts.xlsx (sex, HER2, MSI, EBV, CPS only - never response or survival).
Gene level: TPM and NumReads summed over all transcripts of a gene symbol (transcripts without a symbol keep their Ensembl
         gene ID; alt-haplotype / patch copies, e.g. HLA-A on MHC haplotypes, are summed into their symbol).
Primary expression (decided on QC grounds before the full-cohort analysis; see analysis/results/E22_README.md):
         TPM re-normalised to 1e6 over NUCLEAR protein-coding + IG/TR gene segments ("nuclear coding"), log2(TPM + 1) once.
         Library batch = adapter chemistry measured from the reads (TruSeq-type vs Nextera-type).
Outputs (gene symbol x sample E##):
         PRJEB60680_baseline_gene_TPM.tsv.gz              all quantified transcripts (as quantified)
         PRJEB60680_baseline_gene_TPM_nuclear_coding.tsv.gz   primary TPM
         PRJEB60680_baseline_log2TPM1.tsv.gz              primary log2(TPM + 1)  (nuclear coding)
         PRJEB60680_baseline_log2TPM1_alltranscripts.tsv.gz
         PRJEB60680_baseline_gene_counts.tsv.gz           NumReads summed by symbol (all transcripts)
         PRJEB60680_baseline_gene_TPM_untrimmed[_nuclear_coding].tsv.gz, PRJEB60680_E36_single_mate_gene_TPM[_nuclear_coding].tsv.gz
         PRJEB60680_baseline_qc.tsv, PRJEB60680_qc_identity_checks.tsv, PRJEB60680_qc_summary.txt, PRJEB60680_qc.png
"""
import os, gzip, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

B = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/validation/bulk_PRJEB60680_Samsung_NCT04249739"
QD, FQ = f"{B}/quant", f"{B}/fastq"
QOUT = os.environ.get("E22_QC_OUT", QD)           # override only for test runs
lines = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); lines.append(s)

tx = pd.read_csv(f"{QD}/tx2gene_Ensembl110_cDNA.tsv.gz", sep="\t")
t2g = dict(zip(tx.transcript_id, tx.gene_name))
gtab = tx.groupby("gene_name").agg(biotypes=("gene_biotype", lambda x: set(x)), chrom=("chrom", lambda x: set(x)))
CODING = {g for g, r in gtab.iterrows()
          if (("protein_coding" in r.biotypes) or any(b.startswith(("IG_", "TR_")) and b.endswith("_gene") for b in r.biotypes))
          and "MT" not in r.chrom}
MITO = {g for g, r in gtab.iterrows() if "MT" in r.chrom}
def nc(T):
    """re-normalise TPM over nuclear protein-coding + IG/TR genes"""
    X = T.loc[[g for g in T.index if g in CODING]]
    return X / X.sum() * 1e6

files = pd.read_csv(f"{FQ}/baseline_files.tsv", sep="\t", header=None, names=["run", "title", "file", "url", "md5", "bytes"])
files["sample"] = files.title.str.extract(r"CO_Pem_(E\d+)_B_TRNA")[0]
inst = pd.read_csv(f"{FQ}/baseline_instrument.tsv", sep="\t", header=None, names=["sample", "instrument", "flowcell"]).set_index("sample")
clin = pd.read_excel(f"{B}/Wright2026_CellRepMed_TableS1_clinical_47pts.xlsx")
clin.columns = ["id", "scr", "sex", "age", "BOR", "fast", "pfs_ev", "os_ev", "pfs_d", "os_d", "HER2", "EBV", "PDL1",
                "cps_B", "cps_F1", "cps_F2", "MSI", "TCGA", "scRNA"]
clin = clin.set_index("id")[["sex", "HER2", "EBV", "MSI", "TCGA", "cps_B"]]     # identity fields only
md5 = pd.read_csv(f"{QD}/md5_check.tsv", sep="\t")

samples = sorted(s for s in files["sample"].unique() if os.path.exists(f"{QD}/salmon/{s}/quant.sf"))
miss = sorted(set(files["sample"]) - set(samples))
say(f"samples with paired quantification: {len(samples)}/{files['sample'].nunique()}; missing {miss}"
    + ("" if miss == ["E36"] else "  (TEST RUN - incomplete)"))
NEX, TRU = "CTGTCTCTTATACACATCT", "AGATCGGAAGAGC"
tpm, cnt, tpm_u, qc = {}, {}, {}, []
for s in samples:
    d = f"{QD}/salmon/{s}"
    q = pd.read_csv(f"{d}/quant.sf", sep="\t"); g = q.Name.map(t2g)
    assert g.notna().all(), f"{s}: transcripts missing from tx2gene"
    tpm[s] = q.TPM.groupby(g).sum(); cnt[s] = q.NumReads.groupby(g).sum()
    m = json.load(open(f"{d}/aux_info/meta_info.json")); lf = json.load(open(f"{d}/lib_format_counts.json"))
    runs = files[files["sample"] == s]
    du = f"{QD}/salmon_untrimmed/{s}"
    qu = pd.read_csv(f"{du}/quant.sf", sep="\t"); tpm_u[s] = qu.TPM.groupby(qu.Name.map(t2g)).sum()
    mu = json.load(open(f"{du}/aux_info/meta_info.json"))
    fp = [json.load(open(f"{QD}/fastp/{r}.json")) for r in sorted(runs.run.unique())]
    raw_pairs = sum(j["summary"]["before_filtering"]["total_reads"] for j in fp) // 2
    # read length and adapter chemistry from the first 200,000 R1 reads (FASTQ still on disk at this stage)
    L, nex, tru = [], 0, 0
    fn = runs[runs.file.str.endswith("_1.fastq.gz")].file.iloc[0]
    with gzip.open(f"{FQ}/{fn}", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                L.append(len(line.strip())); nex += NEX in line; tru += TRU in line
            if i >= 800000:
                break
    mm = md5[md5.file.isin(runs.file)]
    qc.append(dict(sample=s, runs=",".join(sorted(runs.run.unique())), n_runs=runs.run.nunique(), title=runs.title.iloc[0],
                   instrument=inst.instrument.get(s, ""), md5_ok=f"{(mm.status == 'OK').sum()}/{len(runs)}",
                   read_len=int(np.median(L)), min_read_len=int(np.min(L)),
                   nextera_adapter_per_100k=100000 * nex / len(L), truseq_adapter_per_100k=100000 * tru / len(L),
                   raw_pairs=raw_pairs, pairs_after_trim=sum(j["summary"]["after_filtering"]["total_reads"] for j in fp) // 2,
                   adapter_trimmed_read_frac=sum(j["adapter_cutting"]["adapter_trimmed_reads"] for j in fp if "adapter_cutting" in j) / (2 * raw_pairs),
                   insert_peak=fp[0]["insert_size"]["peak"], dup_rate=fp[0]["duplication"]["rate"],
                   q30_rate=fp[0]["summary"]["before_filtering"]["q30_rate"],
                   fragments=m["num_processed"], mapped=m["num_mapped"], mapping_rate=m["percent_mapped"],
                   mapped_of_raw_pairs=100 * m["num_mapped"] / raw_pairs, mapping_rate_untrimmed=mu["percent_mapped"],
                   lib_type=",".join(m["library_types"]), expected_format=lf["expected_format"],
                   compatible_ratio=lf["compatible_fragment_ratio"], frag_len_mean=m["frag_length_mean"],
                   index_hash=m["index_seq_hash"][:12]))
qc = pd.DataFrame(qc).set_index("sample")
qc["lib_batch"] = np.where(qc.nextera_adapter_per_100k > qc.truseq_adapter_per_100k, "Nextera-type", "TruSeq-type")
TPM = pd.DataFrame(tpm).fillna(0); CNT = pd.DataFrame(cnt).fillna(0); TPMU = pd.DataFrame(tpm_u).fillna(0)
TPMNC, TPMUNC = nc(TPM), nc(TPMU)
L2 = np.log2(TPMNC + 1); L2all = np.log2(TPM + 1)
qc["mito_TPM_pct"] = TPM.loc[[g for g in TPM.index if g in MITO]].sum() / 1e4
qc["noncoding_nuclear_TPM_pct"] = 100 - qc.mito_TPM_pct - TPM.loc[TPMNC.index].sum() / 1e4
qc["genes_TPM_ge1"] = (TPMNC >= 1).sum()
say(f"genes (symbols / IDs): {TPM.shape[0]} quantified; {TPMNC.shape[0]} nuclear protein-coding + IG/TR (primary set)")
say(f"TPM column sums: {TPM.sum().min():.0f} - {TPM.sum().max():.0f} (all) ; primary re-normalised to 1e6")

# ---- library batch
say("\n== library batch (adapter chemistry measured in the first 200,000 R1 reads)")
say(pd.crosstab([qc.instrument, qc.read_len], qc.lib_batch).to_string())
for bt, sub in qc.groupby("lib_batch"):
    say(f"[{bt}, n={len(sub)}] Nextera adapter /100k reads median {sub.nextera_adapter_per_100k.median():.0f}, TruSeq {sub.truseq_adapter_per_100k.median():.0f}; "
        f"mito TPM% median {sub.mito_TPM_pct.median():.1f} ({sub.mito_TPM_pct.min():.1f}-{sub.mito_TPM_pct.max():.1f}); "
        f"non-coding nuclear TPM% median {sub.noncoding_nuclear_TPM_pct.median():.1f} (max {sub.noncoding_nuclear_TPM_pct.max():.1f})")

# ---- sample - patient matching
say("\n== sample-patient matching")
say(f"run titles -> E## : {len(samples)} baseline samples with paired data; in clinical table: {sum(s in clin.index for s in samples)}")
say(f"clinical patients without a baseline RNA run: {sorted(set(clin.index) - set(files['sample']))}")
qc["in_clinical"] = [s in clin.index for s in qc.index]
c = clin.reindex(qc.index)
qc["sex_recorded"] = c.sex
qc["RPS4Y1_log2TPM"] = L2.loc["RPS4Y1"]
qc["sex_expr"] = np.where(qc.RPS4Y1_log2TPM > 0.5, "M", "F")      # XIST (lncRNA) is not in the cDNA index
qc["sex_match"] = qc.sex_expr == qc.sex_recorded
say(f"RPS4Y1 log2(TPM+1): recorded M min {qc.RPS4Y1_log2TPM[qc.sex_recorded == 'M'].min():.2f}, recorded F max "
    f"{qc.RPS4Y1_log2TPM[qc.sex_recorded == 'F'].max():.2f}; agreement {qc.sex_match.sum()}/{len(qc)}"
    + ("" if qc.sex_match.all() else f"; MISMATCH: {list(qc.index[~qc.sex_match])}"))
# HER2 / MSI / CPS checks on within-batch ranks (absolute levels differ between chemistries)
Lc = L2.sub(L2.T.groupby(qc.lib_batch).transform("mean").T)            # within-batch centred
idc = []
for gname, lab in [("ERBB2", "HER2"), ("GRB7", "HER2")]:
    v = Lc.loc[gname]; pos = c.HER2 == "positive"
    idc.append(dict(check=f"{gname} (within-batch centred log2TPM), HER2+ vs HER2-", n_pos=int(pos.sum()), med_pos=v[pos].median(),
                    med_neg=v[~pos].median(), p=mannwhitneyu(v[pos], v[~pos]).pvalue,
                    detail="HER2+ rank (1 = highest of %d): " % len(v) + ",".join(f"{s}:{int(v.rank(ascending=False)[s])}" for s in pos[pos].index)))
v = Lc.loc["MLH1"]; msi = c.MSI == "MSI-H"
idc.append(dict(check="MLH1 (within-batch centred), MSI-H vs MSS/MSI-L", n_pos=int(msi.sum()), med_pos=v[msi].median(),
                med_neg=v[~msi].median(), p=mannwhitneyu(v[msi], v[~msi]).pvalue,
                detail="rank from lowest: " + ",".join(f"{s}:{int(v.rank()[s])}" for s in list(msi[msi].index) + list(c.index[c.MSI == 'MSI-L']))))
cps = pd.to_numeric(c.cps_B, errors="coerce"); r, p = spearmanr(cps, Lc.loc["CD274"], nan_policy="omit")
idc.append(dict(check="CD274 (within-batch centred) vs baseline CPS, Spearman", n_pos=int(cps.notna().sum()), med_pos=r, med_neg=np.nan, p=p, detail=""))
ebv = c.EBV == "positive"
idc.append(dict(check="EBV+ (no viral reference in the index; not checkable)", n_pos=int(ebv.sum()), med_pos=np.nan, med_neg=np.nan,
                p=np.nan, detail=",".join(ebv[ebv].index)))
idc = pd.DataFrame(idc)
for _, r in idc.iterrows():
    say(f"  {r.check}: n+={r.n_pos} {r.med_pos:.2f} vs {r.med_neg:.2f}  P={r.p:.3g}  {r.detail}")

# ---- mapping / library
say("\n== mapping and library")
for bt, sub in [("all", qc)] + list(qc.groupby("lib_batch")):
    say(f"[{bt}, n={len(sub)}] mapping % (trimmed): median {sub.mapping_rate.median():.1f} ({sub.mapping_rate.min():.1f}-{sub.mapping_rate.max():.1f}); "
        f"of raw pairs {sub.mapped_of_raw_pairs.median():.1f}; untrimmed {sub.mapping_rate_untrimmed.median():.1f} "
        f"({sub.mapping_rate_untrimmed.min():.1f}-{sub.mapping_rate_untrimmed.max():.1f}); adapter-trimmed reads "
        f"{100 * sub.adapter_trimmed_read_frac.median():.1f}%; insert peak {sub.insert_peak.median():.0f}; fragments "
        f"{sub.fragments.median() / 1e6:.1f}M ({sub.fragments.min() / 1e6:.1f}-{sub.fragments.max() / 1e6:.1f})")
lu = np.log2(TPMUNC.reindex(L2.index).fillna(0) + 1)
rg = [np.corrcoef(L2[s], lu[s])[0, 1] for s in qc.index]; qc["r_trimmed_vs_untrimmed"] = rg
say(f"per-sample Pearson r, trimmed vs untrimmed (primary scale): median {np.median(rg):.4f}, min {np.min(rg):.4f}")
say(f"library types: {qc.lib_type.value_counts().to_dict()}; compatible ratio min {qc.compatible_ratio.min():.3f}")
say(f"read lengths: {qc.read_len.value_counts().to_dict()}; minimum read length before trimming {qc.min_read_len.min()} "
    f"(index k = 31; fastp drops reads < 31 after trimming)")
say(f"md5: {md5.status.value_counts().to_dict()} over {len(md5)} files")
say(f"nuclear coding genes with TPM >= 1: median {qc.genes_TPM_ge1.median():.0f}, range {qc.genes_TPM_ge1.min()}-{qc.genes_TPM_ge1.max()}")

# ---- correlation outliers (within batch) and batch structure
allx = L2[(L2.mean(1) >= 1)]
R = pd.DataFrame(np.corrcoef(allx.T.values), index=L2.columns, columns=L2.columns)
np.fill_diagonal(R.values, np.nan)
same = pd.DataFrame(qc.lib_batch.values[:, None] == qc.lib_batch.values[None, :], index=R.index, columns=R.columns)
qc["median_r_within_batch"] = R.where(same).median(1); qc["median_r_between_batch"] = R.where(~same).median(1)
qc["max_r"] = R.max(1); qc["max_r_partner"] = R.idxmax(1)
qc["corr_outlier"] = False
for bt, sub in qc.groupby("lib_batch"):
    med, mad = sub.median_r_within_batch.median(), (sub.median_r_within_batch - sub.median_r_within_batch.median()).abs().median() * 1.4826
    qc.loc[sub.index, "corr_outlier"] = sub.median_r_within_batch < med - 3 * mad
    say(f"\n[{bt}] median within-batch r {med:.3f} (MAD {mad:.3f}); outliers (< median - 3 MAD): {list(sub.index[sub.median_r_within_batch < med - 3 * mad])}")
say(f"median between-batch r {np.nanmedian(R.where(~same).values):.3f}; highest pair r {np.nanmax(R.values):.3f} {R.stack().idxmax()} "
    f"(duplicate-library flag if > 0.99: {'YES' if np.nanmax(R.values) > 0.99 else 'none'})")
mr_med, mr_mad = qc.mapping_rate.median(), (qc.mapping_rate - qc.mapping_rate.median()).abs().median() * 1.4826
qc["mapping_outlier"] = (qc.mapping_rate < 50) | (qc.mapping_rate < mr_med - 3 * mr_mad)
say(f"mapping-rate outliers (< 50% or < median - 3 MAD): {list(qc.index[qc.mapping_outlier])}")
def pca(M):
    e = M.loc[M.var(1).sort_values(ascending=False).index[:5000]]
    X = e.T.values - e.T.values.mean(0); U, S, Vt = np.linalg.svd(X, full_matrices=False)
    return U[:, :2] * S[:2], S ** 2 / (S ** 2).sum()
nx = qc.lib_batch == "Nextera-type"
for lab, M in [("primary log2 TPM", allx), ("primary, centred within batch", allx.sub(allx.T.groupby(qc.lib_batch).transform("mean").T))]:
    pcs, ve = pca(M)
    if lab.startswith("primary log2"):
        qc["PC1"], qc["PC2"] = pcs[:, 0], pcs[:, 1]; ve_raw = ve
    else:
        qc["PC1_within"], qc["PC2_within"] = pcs[:, 0], pcs[:, 1]; ve_w = ve
    say(f"PCA ({lab}, top-5,000 variable expressed genes): PC1 {100 * ve[0]:.1f}%, PC2 {100 * ve[1]:.1f}%; batch vs PC1 Mann-Whitney "
        f"P = {mannwhitneyu(pcs[nx.values, 0], pcs[~nx.values, 0]).pvalue:.3g}; PC1 vs mito Spearman "
        f"{spearmanr(pcs[:, 0], qc.mito_TPM_pct)[0]:.2f}")
qc["QC_pass"] = qc.in_clinical & qc.sex_match & ~qc.corr_outlier & ~qc.mapping_outlier & qc.md5_ok.str.match(r"(\d+)/\1")
qc["exclusion_reason"] = ""
for s in qc.index[~qc.QC_pass]:
    qc.loc[s, "exclusion_reason"] = ";".join(k for k, bad in [("not in clinical table", not qc.in_clinical[s]),
        ("sex mismatch", not qc.sex_match[s]), ("correlation outlier", qc.corr_outlier[s]), ("mapping outlier", qc.mapping_outlier[s])] if bad)

# ---- E36 (ERR13381252): the deposited mate files come from two different libraries
say("\n== E36 (ERR13381252)")
E36 = {}
for m in ("1", "2"):
    q = pd.read_csv(f"{QD}/E36_mate_check/salmon_R{m}_SE/quant.sf", sep="\t"); E36[f"E36_R{m}only"] = q.TPM.groupby(q.Name.map(t2g)).sum()
    mi = json.load(open(f"{QD}/E36_mate_check/salmon_R{m}_SE/aux_info/meta_info.json"))
    say(f"E36 R{m} alone (single-end, fastp-trimmed): {mi['num_processed']} reads, mapping {mi['percent_mapped']:.1f}%")
E36 = pd.DataFrame(E36).reindex(TPM.index).fillna(0); E36NC = nc(E36); le = np.log2(E36NC + 1)
Lnx = L2.loc[:, nx]                                   # E36 R1/R2 were sequenced on the NovaSeq (Nextera-type) flowcell
AMP = ["ERBB2", "GRB7", "PGAP3", "MIEN1"]
for col in E36.columns:
    pct = {g: int(round(100 * (Lnx.loc[g] < le.loc[g, col]).mean())) for g in AMP}
    say(f"  {col}: RPS4Y1 {le.loc['RPS4Y1', col]:.2f} (male if > 0.5); ERBB2-amplicon percentile vs Nextera-batch samples {pct}; "
        f"nuclear coding genes TPM >= 1: {(E36NC[col] >= 1).sum()} (cohort median {int(qc.genes_TPM_ge1.median())}); "
        f"mito TPM% {E36.loc[[g for g in E36.index if g in MITO], col].sum() / 1e4:.0f}")
    idc = pd.concat([idc, pd.DataFrame([dict(check=f"{col}: RPS4Y1 log2TPM; ERBB2-amplicon percentiles vs Nextera batch", n_pos=1,
                     med_pos=le.loc["RPS4Y1", col], med_neg=np.nan, p=np.nan, detail=str(pct))])], ignore_index=True)
ex = allx.index
say(f"  R1-only vs R2-only Pearson r (expressed genes): {np.corrcoef(le.loc[ex].iloc[:, 0], le.loc[ex].iloc[:, 1])[0, 1]:.3f}; "
    f"each vs Nextera-batch samples: median r {np.median([np.corrcoef(le.loc[ex].iloc[:, 0], L2.loc[ex, s])[0, 1] for s in Lnx.columns]):.3f} / "
    f"{np.median([np.corrcoef(le.loc[ex].iloc[:, 1], L2.loc[ex, s])[0, 1] for s in Lnx.columns]):.3f}")
say("  E36 recorded: male, HER2-positive. Index barcodes differ between the files (R1 TGGTGTTATG+TAGAATTGGA, R2 TGGCCTCTGT+CGCTGTCTCA;"
    " neither occurs in any other RNA or WXS run of PRJEB60680) and no read name is shared, so ENA could not build pairs.")
say("  -> E36 EXCLUDED from the primary analysis; its R2-only quantification (male, ERBB2 amplicon high) is kept for a sensitivity analysis.")
row = pd.DataFrame([dict(runs="ERR13381252", n_runs=1, title="CO_Pem_E36_B_TRNA", instrument="A01192", md5_ok="2/2", in_clinical=True,
                         sex_recorded="M", lib_batch="Nextera-type", QC_pass=False,
                         exclusion_reason="mate files from two different libraries (different index barcodes, no shared read names); "
                                          "paired quantification impossible")], index=["E36"])
qc = pd.concat([qc, row]); qc.index.name = "sample"
say(f"\nQC pass: {int(qc.QC_pass.sum())}/{len(qc)}; excluded: " + "; ".join(f"{s} ({qc.exclusion_reason[s]})" for s in qc.index[~qc.QC_pass.astype(bool)]))
say("batch of QC-pass samples: " + str(qc[qc.QC_pass.astype(bool)].lib_batch.value_counts().to_dict()))

# ---- write
os.chdir(QOUT)
TPM.to_csv("PRJEB60680_baseline_gene_TPM.tsv.gz", sep="\t", float_format="%.4f")
TPMNC.to_csv("PRJEB60680_baseline_gene_TPM_nuclear_coding.tsv.gz", sep="\t", float_format="%.4f")
L2.to_csv("PRJEB60680_baseline_log2TPM1.tsv.gz", sep="\t", float_format="%.5f")
L2all.to_csv("PRJEB60680_baseline_log2TPM1_alltranscripts.tsv.gz", sep="\t", float_format="%.5f")
CNT.to_csv("PRJEB60680_baseline_gene_counts.tsv.gz", sep="\t", float_format="%.3f")
TPMU.to_csv("PRJEB60680_baseline_gene_TPM_untrimmed.tsv.gz", sep="\t", float_format="%.4f")
TPMUNC.to_csv("PRJEB60680_baseline_gene_TPM_untrimmed_nuclear_coding.tsv.gz", sep="\t", float_format="%.4f")
E36.to_csv("PRJEB60680_E36_single_mate_gene_TPM.tsv.gz", sep="\t", float_format="%.4f")
E36NC.to_csv("PRJEB60680_E36_single_mate_gene_TPM_nuclear_coding.tsv.gz", sep="\t", float_format="%.4f")
qc.to_csv("PRJEB60680_baseline_qc.tsv", sep="\t")
idc.to_csv("PRJEB60680_qc_identity_checks.tsv", sep="\t", index=False)
open("PRJEB60680_qc_summary.txt", "w").write("\n".join(lines) + "\n")

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
qp = qc[qc.mapping_rate.notna()]
col = qp.lib_batch.map({"TruSeq-type": "#0F6FC5", "Nextera-type": "#C96F12"})
fig, ax = plt.subplots(1, 4, figsize=(17, 4))
o = qp.sort_values(["lib_batch", "mapping_rate"]).index
ax[0].bar(range(len(o)), qp.loc[o, "mapping_rate"], color=col[o], label="trimmed (primary); blue TruSeq-type, orange Nextera-type")
ax[0].plot(range(len(o)), qp.loc[o, "mapping_rate_untrimmed"], "k_", ms=6, label="untrimmed")
ax[0].set_xticks(range(len(o))); ax[0].set_xticklabels(o, rotation=90, fontsize=5); ax[0].set_ylabel("salmon mapping rate (%)"); ax[0].legend(fontsize=6)
for k, (a, b, t) in enumerate([("PC1", "PC2", f"PCA, primary log2 TPM (PC1 {100 * ve_raw[0]:.0f}%)"),
                               ("PC1_within", "PC2_within", f"PCA, centred within batch (PC1 {100 * ve_w[0]:.0f}%)")]):
    ax[k + 1].scatter(qp[a], qp[b], c=col, s=14)
    for s in qp.index:
        ax[k + 1].annotate(s, (qp[a][s], qp[b][s]), fontsize=5)
    ax[k + 1].set_title(t, fontsize=8); ax[k + 1].set_xlabel(a); ax[k + 1].set_ylabel(b)
o = qp.sort_values("lib_batch").index
im = ax[3].imshow(R.loc[o, o].values, cmap="viridis"); ax[3].set_xticks([]); ax[3].set_yticks(range(len(o)))
ax[3].set_yticklabels(o, fontsize=4); plt.colorbar(im, ax=ax[3], fraction=0.046); ax[3].set_title("Pearson r (sorted by batch)", fontsize=8)
fig.tight_layout(); fig.savefig("PRJEB60680_qc.png", dpi=200)
say("written: matrices, QC table, identity checks, summary, PRJEB60680_qc.png")
