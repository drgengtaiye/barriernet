# v33 Fig. 2 (v31 Fig. 3): panel d shows the first three locked cohorts both ways

Third pre-submission review, point 1: keep the value on which the extension was decided beside the duplicate-free one.
Panel d now has two pooled rows for the first three locked cohorts: "as first analyzed" (with the eight GSE26253
arrays later found to be shared with GSE62254; `analysis/results/E30_dedup_runs/off/E20_pooled.csv`; NetImmune 0.73,
0.53–1.02, P = 0.057) and "without shared patients" (E30 primary; 0.72, 0.53–0.996, P = 0.049). Confidence limits that
would print as 1.00 on the other side of 1 are printed to three decimals. The small arrow under the first-three pool is
removed. Everything else is drawn as in v32 (same overlay of results as `../v32_dedup/`); pixel differences against
the v32 render are confined to panel d.

`make_fig3_v33.py` is the script as run (sys.path and the style module point to a scratch overlay, as in v32);
`make_fig3_v33.diff` is its difference from `../make_fig3.py`, which was not changed.
