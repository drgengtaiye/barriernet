"""Build the lean manuscript (MANUSCRIPT_lean.{md,docx}) and its Supplementary Information
(SUPPLEMENTARY_lean.{md,docx}) from lean_src.md and lean_supp_src.md; each is numbered separately.

Citations in the source are written [@key] or [@key1; @key2]; they are numbered by first
appearance, rendered as superscripts (ranges collapsed), and the reference list is inserted
at {{REFERENCES}}. Entries shared with MANUSCRIPT_GM.md are taken verbatim from its list.
"""
import re, subprocess, pathlib
from refs import R, HERE, OUT

def fmt(ns):
    ns = sorted(set(ns)); out = []; i = 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1: j += 1
        out.append(str(ns[i]) if j == i else (f"{ns[i]},{ns[j]}" if j == i + 1 else f"{ns[i]}–{ns[j]}"))
        i = j + 1
    return "^" + ",".join(out) + "^"


def build(src_name, out_name, label):
    src = (HERE / src_name).read_text()
    body, tail = src.split("{{REFERENCES}}")
    cite = re.compile(r"\[(@[^\]]+)\]")
    order = []
    for m in cite.finditer(body):
        for k in re.findall(r"@(\w+)", m.group(1)):
            if k not in order: order.append(k)
    missing = [k for k in order if k not in R]; assert not missing, missing
    num = {k: i + 1 for i, k in enumerate(order)}
    body = cite.sub(lambda m: fmt([num[k] for k in re.findall(r"@(\w+)", m.group(1))]), body)
    out = body + "\n\n".join(f"{num[k]}\\. {R[k]}" for k in order) + tail
    (OUT / f"{out_name}.md").write_text(out)
    subprocess.run(["pandoc", f"{out_name}.md", "-o", f"{out_name}.docx"], check=True, cwd=OUT)
    return out, order


main, order_main = build("lean_src.md", "MANUSCRIPT_lean", "main")
supp, order_supp = build("lean_supp_src.md", "SUPPLEMENTARY_lean", "supp")
abstract = main.split("## Abstract")[1].split("**Keywords")[0]
def words(t): return len(t.split())
intro_to_disc = main.split("## Introduction")[1].split("## Declarations")[0]
tables = main.split("## Tables")[1].split("## References")[0]
legends = main.split("## Figures")[1]
print(f"main: {len(order_main)} references | abstract {words(abstract)} words (incl. 4 labels) | "
      f"Introduction-Conclusion {words(intro_to_disc)} words (text only) | tables {words(tables)} | legends {words(legends)}")
print(f"supplement: {len(order_supp)} references | {words(supp)} words")
print("keys in R not cited in main:", [k for k in R if k not in order_main])
