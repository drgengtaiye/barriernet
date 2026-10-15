"""PRJEB60680 baseline bulk RNA-seq: md5-verify each downloaded FASTQ pair, trim, then salmon-quantify it.

Runs on the compute node while login1 downloads (fastq/dl_all.sh). Polls fastq/ until every baseline run listed
in fastq/baseline_files.tsv has both mates at the expected byte size, checks md5 against ENA (fastq_md5, or
submitted_md5 for ERR13381252 whose ENA fastq is not generated), then per SAMPLE (all runs of one patient's
baseline biopsy passed together, i.e. merged; here one run each):
  1. untrimmed (sensitivity; v1 of this script):  salmon quant -l A -1 R1 -2 R2 --validateMappings -> salmon_untrimmed/<E##>
  2. primary: fastp 0.24.0 (PE overlap trimming + --adapter_fasta adapters.fa [TruSeq R1/R2, Nextera/Tn5 mosaic end]
     --trim_poly_g --length_required 31, per run; sequence auto-detection NOT used because on a 2x101 sample it
     reported transcript sequence as "adapter") -> salmon quant
     -l A --validateMappings on the trimmed pairs -> salmon/<E##>; trimmed FASTQ deleted after quantification.
     Reason (found 2026-10-05 after the first 21 samples): the 2x151 NovaSeq libraries carry Nextera/Tn5 adapter
     read-through (CTGTCTCTTATACACATCT in ~41% of R1, poly-G tails ~15%; insert peak ~120 bp), so selective
     alignment discarded ~20% of fragments (mapping 47-56% vs 66-77% in the 2x101 batch). Trimming is applied
     identically to every sample.
Index: Ensembl 110 GRCh38 cDNA (k = 31), salmon 2.0.0. 3 samples at a time x 32 threads ~ nproc (101).
"""
import os, time, hashlib, subprocess, shutil, concurrent.futures as cf
import pandas as pd

B = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/validation/bulk_PRJEB60680_Samsung_NCT04249739"
FQ, QD = f"{B}/fastq", f"{B}/quant"
SALMON = "/groups/adv2105_gp/yichen/Yi/com/gastric_xiong/env/salmon/bin/salmon"
FASTP = "/groups/adv2105_gp/yichen/tools/fastp"
IDX = "/groups/adv2105_gp/yichen/Yi/com/gastric_xiong/ref/salmon_idx"
THREADS, PAR = 32, 3
for d in ("salmon", "salmon_untrimmed", "fastp", "trim_tmp"):
    os.makedirs(f"{QD}/{d}", exist_ok=True)

f = pd.read_csv(f"{FQ}/baseline_files.tsv", sep="\t", header=None,
                names=["run", "title", "file", "url", "md5", "bytes"])
f["sample"] = f.title.str.extract(r"CO_Pem_(E\d+)_B_TRNA")[0]
f["mate"] = f.file.str.extract(r"_([12])\.fastq\.gz")[0]
samples = sorted(f["sample"].unique())
MD5LOG = f"{QD}/md5_check.tsv"
done_md5 = {}
if os.path.exists(MD5LOG):
    for line in open(MD5LOG):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0] != "file":
            done_md5[p[0]] = p[3]
else:
    open(MD5LOG, "w").write("file\texpected\tobserved\tstatus\n")


def log(msg):
    with open(f"{QD}/quant_watch.log", "a") as fh:
        fh.write(time.strftime("%F %T ") + msg + "\n")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def ready(sample):
    for _, r in f[f["sample"] == sample].iterrows():
        p = f"{FQ}/{r.file}"
        if not os.path.exists(p) or os.path.getsize(p) != r.bytes:
            return False
    return True


def salmon(r1, r2, out, logf):
    cmd = [SALMON, "quant", "-i", IDX, "-l", "A", "-1", *r1, "-2", *r2, "--validateMappings", "-p", str(THREADS), "-o", out]
    with open(logf, "w") as lf:
        rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
    return rc == 0 and os.path.exists(f"{out}/quant.sf")


def process(sample):
    sub = f[f["sample"] == sample]
    for _, r in sub.iterrows():
        if done_md5.get(r.file) == "OK":
            continue
        obs = md5(f"{FQ}/{r.file}")
        st = "OK" if obs == r.md5 else "MISMATCH"
        with open(MD5LOG, "a") as fh:
            fh.write(f"{r.file}\t{r.md5}\t{obs}\t{st}\n")
        done_md5[r.file] = st
        if st != "OK":
            log(f"{sample} md5 MISMATCH {r.file}")
            return sample, "md5_fail"
    runs = sorted(sub.run.unique())
    r1 = [f"{FQ}/{run}_1.fastq.gz" for run in runs]; r2 = [f"{FQ}/{run}_2.fastq.gz" for run in runs]
    t0 = time.time()
    if not os.path.exists(f"{QD}/salmon_untrimmed/{sample}/quant.sf"):
        ok = salmon(r1, r2, f"{QD}/salmon_untrimmed/{sample}", f"{QD}/salmon_untrimmed/{sample}.log")
        log(f"{sample} untrimmed salmon {'OK' if ok else 'FAIL'} {time.time() - t0:.0f}s")
        if not ok:
            return sample, "salmon_fail"
    if os.path.exists(f"{QD}/salmon/{sample}/quant.sf"):
        return sample, "skip"
    t1, t2 = [], []
    for run in runs:
        o1, o2 = f"{QD}/trim_tmp/{run}_1.fq.gz", f"{QD}/trim_tmp/{run}_2.fq.gz"
        cmd = [FASTP, "-i", f"{FQ}/{run}_1.fastq.gz", "-I", f"{FQ}/{run}_2.fastq.gz", "-o", o1, "-O", o2,
               "--adapter_fasta", f"{QD}/adapters.fa", "--trim_poly_g", "--length_required", "31", "-w", "16", "-z", "1",
               "-j", f"{QD}/fastp/{run}.json", "-h", f"{QD}/fastp/{run}.html"]
        with open(f"{QD}/fastp/{run}.log", "w") as lf:
            rc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            log(f"{sample} fastp FAIL {run}")
            return sample, "fastp_fail"
        t1.append(o1); t2.append(o2)
    ok = salmon(t1, t2, f"{QD}/salmon/{sample}", f"{QD}/salmon/{sample}.log")
    for x in t1 + t2:
        os.remove(x)
    log(f"{sample} trimmed salmon {'OK' if ok else 'FAIL'} {time.time() - t0:.0f}s runs={','.join(runs)}")
    return sample, "ok" if ok else "salmon_fail"


status = {s: ("done" if os.path.exists(f"{QD}/salmon/{s}/quant.sf") else "pending") for s in samples}
log(f"start v2 (trimmed primary): {sum(v == 'done' for v in status.values())}/{len(samples)} already quantified")
running = {}
with cf.ThreadPoolExecutor(max_workers=PAR) as ex:
    while True:
        for fut in [x for x in running if x.done()]:
            s, st = fut.result(); status[s] = st; del running[fut]
        for s in [s for s in samples if status[s] == "pending"]:
            if len(running) >= PAR:
                break
            if ready(s):
                status[s] = "running"; running[ex.submit(process, s)] = s
        if not running and all(status[s] not in ("pending", "running") for s in samples):
            break
        if not running and os.path.exists(f"{QD}/STOP"):
            break
        time.sleep(60)
log("finished: " + ", ".join(f"{k}={v}" for k, v in sorted(status.items()) if v not in ("ok", "done", "skip")))
log("all statuses: " + str(pd.Series(status).value_counts().to_dict()))
