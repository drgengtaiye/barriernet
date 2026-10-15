"""Assemble the submission package in final/ from the lean build (run build_lean.py first).

final/
  BarrierNet_Manuscript.docx                  title page, main text, references, tables, figures with legends
  BarrierNet_Supplementary_Information.docx   Supplementary Methods, Tables S2-S4/S8/S9, Figures S1-S20
  BarrierNet_Supplementary_Tables.xlsx        workbook (Supplementary Tables S1, S5-S7, S10-S15)
  figures/Figure1-5.{pdf,png}, figures/supplementary/FigureS1-S20.{pdf,png}

Layout: US Letter, 1-inch margins, Times New Roman 12 pt, double-spaced main text with continuous line
numbers, page numbers, References before Tables, table notes below each table in 9 pt, and each table and
figure on its own page. The text itself is taken unchanged from MANUSCRIPT_lean.md / SUPPLEMENTARY_lean.md.
"""
import re, shutil, subprocess, pathlib
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Inches, RGBColor

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
import sys
VER = sys.argv[1] if len(sys.argv) > 1 else "lean"          # "lean" -> final/, "v18" -> final_v18/
FINAL = ROOT / ("final" if VER == "lean" else f"final_{VER}")
FIGDIR = ROOT / "figures" / VER
WORKBOOK = FIGDIR / f"Supplementary_Tables_{VER}.xlsx"
FONT = "Times New Roman"
RUNNING_TITLE = "barriernet: prognostic versus predictive information"


def sections(md):
    """Split at level-2 headings; returns (preamble, [(heading, text), ...])."""
    parts = re.split(r"(?m)^(## .+)$", md)
    return parts[0], [(parts[i][3:].strip(), parts[i + 1]) for i in range(1, len(parts), 2)]


def words(t):
    t = re.sub(r"\^[0-9,–]+\^", "", t)
    t = "\n".join(l for l in t.splitlines() if not l.startswith("#"))
    return len(re.sub(r"[*_`]", "", t).split())


def notes_below(tables_md):
    """Within each table block (caption, bullet notes, pipe table), move the notes below the table."""
    out = []
    for blk in re.split(r"(?m)^(?=\*\*(?:Supplementary )?Table S?\d+\.)", tables_md):
        lines = blk.strip("\n").splitlines()
        if not lines or not lines[0].startswith("**"):
            out.append(blk.strip("\n")); continue
        cap = lines[0]
        tab = [l for l in lines[1:] if l.startswith("|")]
        notes = [l for l in lines[1:] if l.strip() and not l.startswith("|")]
        out.append("\n".join([cap, ""] + tab + [""] + notes))
    return "\n\n".join(x for x in out if x.strip()) + "\n\n"


# ---------------------------------------------------------------- markdown restructuring
def main_md():
    md = (ROOT / f"MANUSCRIPT_{VER}.md").read_text()
    pre, secs = sections(md)
    S = dict(secs)
    title = pre.strip().splitlines()[0]
    assert title.startswith("# ")
    abstract = S["Abstract"].split("**Keywords")[0]
    body = "".join(S[h] for h in ["Introduction", "Materials and Methods", "Results", "Discussion", "Conclusion"])
    n_ref = len(re.findall(r"(?m)^\d+\\\. ", S["References"]))
    sup = (ROOT / f"SUPPLEMENTARY_{VER}.md").read_text()
    n_sfig = len(set(re.findall(r"\*\*Supplementary Fig\. S(\d+)\.", sup)))
    n_stab = len(set(re.findall(r"Supplementary Table S(\d+)\.", sup)))
    n_fig = len(re.findall(r"(?m)^!\[\]", S["Figures"]))
    n_tab = len(re.findall(r"(?m)^\*\*Table \d+\.", S["Tables"]))
    counts = (f"abstract {words(abstract) - 4} words (excluding the four section labels); main text (Introduction to "
              f"Conclusion) {words(body):,} words; {n_fig} figures; {n_tab} tables; {n_ref} references. "
              f"Supplementary Information: {n_sfig} figures and {n_stab} tables.")
    tp = f"""{title}

**Running title:** {RUNNING_TITLE}

**Authors:** [to be completed]

**Affiliations:** [to be completed]

**Corresponding author:** [name, postal address and e-mail — to be completed]

**Word counts:** {counts}
"""
    order = ["Abstract", "Introduction", "Materials and Methods", "Results", "Discussion", "Conclusion",
             "Declarations", "References", "Tables", "Figures"]
    assert set(order) == set(S), set(S) ^ set(order)
    out = tp
    for h in order:
        txt = S[h]
        if h == "Tables":
            txt = notes_below(txt)
        out += f"\n\n## {h}\n{txt}"
    out = re.sub(r"(?m)^---[ \t]*\n", "", out)          # rules; sections start on new pages instead
    return out, counts


def supp_md():
    md = (ROOT / f"SUPPLEMENTARY_{VER}.md").read_text()
    pre, secs = sections(md)
    out = pre
    for h, txt in secs:
        if h == "Supplementary Tables":
            i = txt.index("### Supplementary Tables in the workbook")
            txt = notes_below(txt[:i]) + txt[i:]
        out += f"\n\n## {h}\n{txt}"
    return re.sub(r"(?m)^---[ \t]*\n", "", out)


# ---------------------------------------------------------------- docx formatting
def set_font(style, size=None, bold=None, italic=None, color=True):
    f = style.font
    f.name = FONT
    rpr = style.element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts"); rpr.insert(0, rf)
    for a in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        rf.set(qn(a), FONT)
    for a in ("w:asciiTheme", "w:hAnsiTheme", "w:cstheme", "w:eastAsiaTheme"):
        if rf.get(qn(a)) is not None:
            del rf.attrib[qn(a)]
    if size: f.size = Pt(size)
    if bold is not None: f.bold = bold
    if italic is not None: f.italic = italic
    if color: f.color.rgb = RGBColor(0, 0, 0)


def para_fmt(style, spacing, before=0, after=6):
    pf = style.paragraph_format
    pf.line_spacing_rule = {2: WD_LINE_SPACING.DOUBLE, 1: WD_LINE_SPACING.SINGLE}.get(spacing, WD_LINE_SPACING.MULTIPLE)
    if spacing not in (1, 2): pf.line_spacing = spacing
    pf.space_before = Pt(before); pf.space_after = Pt(after)


def page_break_before(p):
    p.paragraph_format.page_break_before = True


def border(el_pr, tag, edges):
    b = OxmlElement(tag)
    for e, sz in edges.items():
        x = OxmlElement(f"w:{e}")
        x.set(qn("w:val"), "single" if sz else "nil"); x.set(qn("w:sz"), str(sz or 0))
        x.set(qn("w:space"), "0"); x.set(qn("w:color"), "000000")
        b.append(x)
    old = el_pr.find(qn(tag))
    if old is not None: el_pr.remove(old)
    el_pr.append(b)


def format_tables(doc):
    for t in doc.tables:
        tpr = t._tbl.tblPr
        w = tpr.find(qn("w:tblW"))
        if w is None:
            w = OxmlElement("w:tblW"); tpr.append(w)
        w.set(qn("w:type"), "pct"); w.set(qn("w:w"), "5000")
        border(tpr, "w:tblBorders", {"top": 8, "left": 0, "bottom": 8, "right": 0, "insideH": 0, "insideV": 0})
        # column widths in proportion to content (longest cell, capped), on a fixed layout over 6.5 in
        ncol = len(t.columns)
        # group-header rows (text only in the first cell) span the whole row and do not set widths
        is_group = [ri > 0 and row.cells[0].text.strip() and not any(c.text.strip() for c in row.cells[1:ncol])
                    for ri, row in enumerate(t.rows)]
        need = [0] * ncol
        for ri, row in enumerate(t.rows):
            if is_group[ri]: continue
            for j, c in enumerate(row.cells[:ncol]):
                need[j] = max(need[j], min(len(c.text.strip()), 30))
        need = [n + 4 for n in need]
        widths = [9360 * n / sum(need) for n in need]
        longest = [0] * ncol                                   # keep the longest word (or a short cell) on one line
        for ri, row in enumerate(t.rows):
            if is_group[ri]: continue
            for j, c in enumerate(row.cells[:ncol]):
                txt = c.text.strip()
                longest[j] = max([longest[j], len(txt) if len(txt) <= 10 else 0] + [len(w) for w in txt.split()])
        floor = [min(95 * L + 200, 3000) for L in longest]
        total = 9360
        if sum(floor) >= total:                              # not enough room: shrink the floors proportionally
            widths = [total * f / sum(floor) for f in floor]
        else:                                                # water-filling: fix columns at their floor, share the rest by need
            fixed = {}
            while True:
                free = [j for j in range(ncol) if j not in fixed]
                room = total - sum(fixed.values())
                widths = [fixed.get(j, room * need[j] / sum(need[k] for k in free)) for j in range(ncol)]
                low = [j for j in free if widths[j] < floor[j]]
                if not low: break
                for j in low: fixed[j] = floor[j]
        widths = [int(w) for w in widths]
        lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed"); tpr.append(lay)
        grid = t._tbl.tblGrid
        for gc, wd in zip(grid.findall(qn("w:gridCol")), widths):
            gc.set(qn("w:w"), str(wd))
        for ri, row in enumerate(t.rows):
            for c, wd in zip(row.cells, widths):
                tcpr = c._tc.get_or_add_tcPr()
                tcw = tcpr.find(qn("w:tcW"))
                if tcw is None:
                    tcw = OxmlElement("w:tcW"); tcpr.insert(0, tcw)
                tcw.set(qn("w:w"), str(wd)); tcw.set(qn("w:type"), "dxa")
        for ri, row in enumerate(t.rows):
            if is_group[ri]:
                mc = row.cells[0].merge(row.cells[ncol - 1])
                tcw = mc._tc.get_or_add_tcPr().find(qn("w:tcW"))
                tcw.set(qn("w:w"), str(sum(widths))); tcw.set(qn("w:type"), "dxa")
        for ri, row in enumerate(t.rows):
            trpr = row._tr.get_or_add_trPr()
            if ri == 0:
                h = OxmlElement("w:tblHeader"); trpr.append(h)
            cant = OxmlElement("w:cantSplit"); trpr.append(cant)
            for c in row.cells:
                if ri == 0:
                    border(c._tc.get_or_add_tcPr(), "w:tcBorders", {"bottom": 6})
                for p in c.paragraphs:
                    pf = p.paragraph_format
                    pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
                    pf.space_before = Pt(1); pf.space_after = Pt(1)
                    for r in p.runs:
                        r.font.size = Pt(8.5 if len(row.cells) >= 7 else 9)
                        if ri == 0: r.font.bold = True


def add_page_numbers(section):
    p = section.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run()
    for kind, text in (("begin", None), (None, "PAGE"), ("end", None)):
        if kind:
            e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), kind)
        else:
            e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = text
        r._r.append(e)
    r.font.name = FONT; r.font.size = Pt(10)
    ppr = p._p.get_or_add_pPr()
    ppr.append(OxmlElement("w:suppressLineNumbers"))


def finish(docx_path, body_spacing, line_numbers, break_before):
    doc = Document(docx_path)
    st = doc.styles
    for name in ["Normal", "Body Text", "First Paragraph", "Compact", "Block Text", "Bibliography"]:
        if name in [s.name for s in st]:
            set_font(st[name], 12 if body_spacing == 2 else 11)
            para_fmt(st[name], body_spacing, after=0 if body_spacing == 2 else 6)
    for name, size in [("Heading 1", 16), ("Heading 2", 14), ("Heading 3", 12), ("Heading 4", 12)]:
        set_font(st[name], size, bold=True, italic=(name == "Heading 3" and False))
        para_fmt(st[name], 1, before=12, after=6)
        st[name].paragraph_format.keep_with_next = True
    set_font(st["Hyperlink"], color=False)
    for s in ("Captioned Figure", "Figure"):
        if s in [x.name for x in st]:
            st[s].paragraph_format.keep_with_next = True
            para_fmt(st[s], 1, before=0, after=6)
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Inches(8.5), Inches(11)
    for m in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, m, Inches(1))
    if line_numbers:
        ln = OxmlElement("w:lnNumType")
        ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "360")
        pgsz = sec._sectPr.find(qn("w:pgSz"))
        pgsz.addnext(ln) if pgsz is not None else sec._sectPr.append(ln)
    add_page_numbers(sec)
    format_tables(doc)

    # page breaks, figure sizes, captions
    max_h, max_w = Inches(6.4), Inches(6.5)
    first_after = {}
    paras = doc.paragraphs
    for i, p in enumerate(paras):
        sty = p.style.name
        txt = p.text.strip()
        if sty == "Heading 2" and txt in break_before:
            page_break_before(p); first_after[txt] = i
        in_list = p._p.pPr is not None and p._p.pPr.numPr is not None
        if re.match(r"(Supplementary )?Table S?\d+\.", txt) and not in_list:
            p.paragraph_format.keep_with_next = True
            prev = paras[i - 1]
            if prev.style.name != "Heading 2":
                page_break_before(p)
        shapes = p._p.findall(".//" + qn("wp:inline"))
        if shapes:
            prev = paras[i - 1]
            if prev.style.name != "Heading 2":
                page_break_before(p)
            p.paragraph_format.keep_with_next = True
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            ext = shapes[0].find(qn("wp:extent"))
            cx, cy = int(ext.get("cx")), int(ext.get("cy"))
            s = min(1.0, max_h / cy, max_w / cx)
            ncx, ncy = int(cx * s), int(cy * s)
            ext.set("cx", str(ncx)); ext.set("cy", str(ncy))
            for a in shapes[0].iter(qn("a:ext")):
                a.set("cx", str(ncx)); a.set("cy", str(ncy))
        # figure legends and table notes: smaller, single-spaced
        if re.match(r"(Supplementary )?(Figure|Fig\.) S?\d+\.", txt):
            for r in p.runs: r.font.size = Pt(10)
            p.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
    # table notes (paragraphs between a table and the next caption/heading) in 9 pt
    body = doc.element.body
    after_table = False
    for el in body.iterchildren():
        if el.tag == qn("w:tbl"):
            after_table = True; continue
        if el.tag != qn("w:p"):
            continue
        ps = el.find(qn("w:pPr"))
        sty = ps.find(qn("w:pStyle")).get(qn("w:val")) if ps is not None and ps.find(qn("w:pStyle")) is not None else ""
        text = "".join(t.text or "" for t in el.iter(qn("w:t")))
        if sty.startswith("Heading") or re.match(r"(Supplementary )?Table S?\d+\.", text):
            after_table = False; continue
        if after_table:
            for rr in el.iter(qn("w:r")):
                rpr = rr.find(qn("w:rPr"))
                if rpr is None:
                    rpr = OxmlElement("w:rPr"); rr.insert(0, rpr)
                sz = OxmlElement("w:sz"); sz.set(qn("w:val"), "18"); rpr.append(sz)
            if ps is None:
                ps = OxmlElement("w:pPr"); el.insert(0, ps)
            sp = OxmlElement("w:spacing"); sp.set(qn("w:line"), "240"); sp.set(qn("w:lineRule"), "auto")
            sp.set(qn("w:after"), "40"); ps.append(sp)
    # title on the title page: plain bold, left-aligned
    doc.core_properties.author = ""
    doc.core_properties.title = paras[0].text
    doc.save(docx_path)


def build():
    FINAL.mkdir(exist_ok=True)
    (FINAL / "figures" / "supplementary").mkdir(parents=True, exist_ok=True)
    m, counts = main_md()
    s = supp_md()
    (FINAL / "BarrierNet_Manuscript.md").write_text(m)
    (FINAL / "BarrierNet_Supplementary_Information.md").write_text(s)
    for stem, kw in [("BarrierNet_Manuscript", dict(body_spacing=2, line_numbers=True,
                       break_before={"Abstract", "Introduction", "References", "Tables", "Figures"})),
                     ("BarrierNet_Supplementary_Information", dict(body_spacing=1.15, line_numbers=False,
                       break_before={"Supplementary Methods", "Supplementary Tables", "Supplementary Figures",
                                     "Supplementary references"}))]:
        out = FINAL / f"{stem}.docx"
        subprocess.run(["pandoc", str(FINAL / f"{stem}.md"), "-f", "markdown", "-o", str(out)], check=True, cwd=ROOT)
        finish(out, **kw)
    shutil.copy(WORKBOOK, FINAL / "BarrierNet_Supplementary_Tables.xlsx")
    n_main = len(list(FIGDIR.glob("Fig[0-9].png"))); n_supp = len(list((FIGDIR / "supp").glob("FigS*.png")))
    for i in range(1, n_main + 1):
        for ext in ("pdf", "png"):
            shutil.copy(FIGDIR / f"Fig{i}.{ext}", FINAL / f"figures/Figure{i}.{ext}")
    for i in range(1, n_supp + 1):
        for ext in ("pdf", "png"):
            shutil.copy(FIGDIR / f"supp/FigS{i}.{ext}", FINAL / f"figures/supplementary/FigureS{i}.{ext}")
    print("word counts:", counts)


if __name__ == "__main__":
    build()
