"""One-time conversion of the lean sources into v18 sources with keyed supplementary references.

In v18_src.md and v18_supp_src.md, supplementary items are cited as {SF:key} and {ST:key};
build_v18.py numbers them by first citation in the main text (then in the supplement) and
collapses consecutive numbers into ranges. Figure files are referenced as ![](FIG:key).
"""
import re, pathlib
HERE = pathlib.Path(__file__).resolve().parent

SF = {1: "sc_atlas", 2: "sc_programmes", 3: "sc_composition", 4: "sc_fractions", 5: "sc_share",
      6: "meta_forest", 7: "fitted_models", 8: "upgrades", 9: "locked", 10: "icb11", 11: "tcga",
      12: "perfect_metabric", 13: "visium", 14: "visium_sens", 15: "barrierindex", 16: "moran",
      17: "transport", 18: "imc", 19: "proteomics", 20: "drugs"}
ST = {1: "headtohead_cohort", 2: "headtohead", 3: "incremental", 4: "fitted", 5: "locked", 6: "modules",
      7: "prjeb_nulls", 8: "prjeb_auc", 9: "prjeb_combo", 10: "icb", 11: "perfect_metabric", 12: "tcga",
      13: "trials", 14: "imc", 15: "gdsc"}


def keyed_list(kind, body):
    """'S1, S2' / 'S19 and S20' / 'S14–S17' -> keyed tokens, keeping separators."""
    M = SF if kind == "SF" else ST
    def rng(m):
        a, b = int(m.group(1)), int(m.group(2))
        return ", ".join(f"S{i}" for i in range(a, b + 1))
    body = re.sub(r"S(\d+)–S(\d+)", rng, body)
    return re.sub(r"S(\d+)", lambda m: "{" + kind + ":" + M[int(m.group(1))] + "}", body)


def convert(t):
    sep = r"(?:S\d+(?:–S\d+)?)(?:(?:, | and )S\d+(?:–S\d+)?)*"
    t = re.sub(r"(Supplementary Figs?\.? )(" + sep + ")", lambda m: m.group(1) + keyed_list("SF", m.group(2)), t)
    t = re.sub(r"(Supplementary Tables? )(" + sep + ")", lambda m: m.group(1) + keyed_list("ST", m.group(2)), t)
    t = re.sub(r"!\[\]\(figures/lean/supp/FigS(\d+)\.png\)\{width=100%\}", lambda m: f"![](FIG:{SF[int(m.group(1))]})", t)
    t = re.sub(r"!\[\]\(figures/lean/Fig(\d)\.png\)\{width=100%\}", lambda m: f"![](FIG:main{m.group(1)})", t)
    assert not re.search(r"Supplementary (Figs?\.?|Tables?) S\d", t), re.findall(r"Supplementary (?:Figs?\.?|Tables?) S\d+.{0,20}", t)
    return t


if __name__ == "__main__":
    for a, b in [("lean_src.md", "v18_src.md"), ("lean_supp_src.md", "v18_supp_src.md")]:
        out = HERE / b
        assert not out.exists(), f"{b} exists; this conversion runs once"
        out.write_text(convert((HERE / a).read_text()))
        print("wrote", b)
