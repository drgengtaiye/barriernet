# v32 re-renders without the eight GSE26253 arrays shared with GSE62254 (E30 Q1)

Fig3 (v32 Fig. 2), Fig8 (v32 Fig. 6) and FigS_locked (Supplementary Fig. S10) were drawn by the unchanged scripts
`make_fig3.py`, `make_fig8.py` and `supp/make_supp_gastric.py locked`, run against an overlay of
`analysis/results/` in which the E15, E20 and E29 files are replaced by `analysis/results/E30_dedup_runs/primary/`
(and `figures/cns/data/E20_null_draws.csv` by its duplicate-free version), with `cns_style.RES` and `OUT` pointed at
the overlay (`cns_style_overlay.py`). The only script edit was the hard-coded "2,152 patients" label in make_fig8.py
(2,144), which is not drawn. The v31 renders in figures/cns/ are unchanged. Pixel differences against v31 are
confined to the duplicate-affected rows and labels.
