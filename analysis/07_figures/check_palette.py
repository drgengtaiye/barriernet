"""Python port of the dataviz validate_palette.js core checks (node is not installed here):
OKLab ΔE×100 under Machado-2009 protan/deutan simulation (target ≥ 8, floor 6), unsimulated
ΔE (floor 15), OKLCH chroma floor 0.10, WCAG contrast vs the light surface (≥ 3)."""
import itertools, math, sys
M = {"protan": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
     "deutan": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
     "tritan": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]]}
s2l = lambda c: c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
lin = lambda h: [s2l(int(h.lstrip("#")[i:i + 2], 16) / 255) for i in (0, 2, 4)]
def oklab(rgb):
    r, g, b = rgb
    l, m, s = (math.copysign(abs(x) ** (1 / 3), x) for x in (0.4122214708*r+0.5363325363*g+0.0514459929*b,
               0.2119034982*r+0.6806995451*g+0.1073969566*b, 0.0883024619*r+0.2817188376*g+0.6299787005*b))
    return (0.2104542553*l+0.7936177850*m-0.0040720468*s, 1.9779984951*l-2.4285922050*m+0.4505937099*s,
            0.0259040371*l+0.7827717662*m-0.8086757660*s)
def sim(h, k):
    r = lin(h); return [max(0, min(1, sum(M[k][i][j] * r[j] for j in range(3)))) for i in range(3)]
def dE(a, b, k=None):
    return 100 * math.dist(oklab(sim(a, k) if k else lin(a)), oklab(sim(b, k) if k else lin(b)))
lum = lambda h: sum(w * c for w, c in zip((0.2126, 0.7152, 0.0722), lin(h)))
def contrast(a, b):
    x, y = sorted((lum(a), lum(b)), reverse=True); return (x + 0.05) / (y + 0.05)
pal = sys.argv[1].split(","); ok = True
for h in pal:
    L, a, b = oklab(lin(h)); C = math.hypot(a, b); cr = contrast(h, "#ffffff")
    flag = "" if C >= 0.10 and cr >= 3 else "  <-- chroma/contrast"
    print(f"{h}: L={L:.2f} C={C:.3f} contrast={cr:.2f}{flag}")
for a, b in itertools.combinations(pal, 2):
    cvd = min(dE(a, b, "protan"), dE(a, b, "deutan")); n = dE(a, b)
    verdict = "PASS" if cvd >= 8 and n >= 15 else ("FLOOR (needs 2nd encoding)" if cvd >= 6 and n >= 15 else "FAIL")
    ok &= verdict != "FAIL"
    print(f"{a}-{b}: CVD {cvd:5.1f}  normal {n:5.1f}  tritan {dE(a, b, 'tritan'):5.1f}  {verdict}")
print("OVERALL", "PASS" if ok else "FAIL")
