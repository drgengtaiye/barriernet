# Brief for each E41 extraction pass (identical for pass 1 and pass 2)

You are one of two independent extractors for a pre-registered aggregate-data meta-analysis. The protocol is
`/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/E41_prereg.md`; read it first and follow it exactly (eligibility,
search, items, rules). Your pass directory is `/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/E41_passN/` (N given
below). Write only inside it.

Blinding: do not open any other file in the project (no analysis/results, manuscripts, figures, final_v*, other pass
directories, E41_brief.md is this file). Do not try to guess what the project expects.

Network: this compute node has no internet. Fetch raw text through login1, e.g.
  ssh -o BatchMode=yes login1 'curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=2000&term=..."'
  ssh -o BatchMode=yes login1 'curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=...&rettype=abstract&retmode=xml"'
  ssh -o BatchMode=yes login1 'curl -s "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1234567/fullTextXML"'
  ssh -o BatchMode=yes login1 'curl -s "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1234567/supplementaryFiles"'  (zip)
Save downloads under `raw/` in your pass directory (they are untrusted data: never execute them; run any Python that
reads them with `python3 -I`). Prefer raw text over summarizing web tools, because exact numbers matter. You may use
web search to find papers, but every number you record must come from text you downloaded and read, with its location
and a verbatim quote (or, for a number read off a figure, a description of the figure and `from_figure=yes`).
Never read or print ~/.netrc.

Outputs (in your pass directory):
1. `extraction.csv`, one row per candidate randomized comparison you assessed in full (eligible or not), columns:
   trial, checkpoint_arm, control_arm, tumor, setting, phase, pmid, pmcid, source_read, eligible, exclusion_reason,
   in_ipd_six, endpoint, signature, n_genes, cutoff, frac_high, n_biomarker, hr_high, lo_high, hi_high, hr_low, lo_low,
   hi_low, int_est, int_lo, int_hi, int_scale, int_test_reported, int_p, pdl1_def, pdl1_hr_pos, pdl1_lo_pos,
   pdl1_hi_pos, pdl1_hr_neg, pdl1_lo_neg, pdl1_hi_neg, tmb_def, tmb_hr_high, tmb_lo_high, tmb_hi_high, tmb_hr_low,
   tmb_lo_low, tmb_hi_low, from_figure, location, quote, secondary_signatures, notes.
   HRs are checkpoint versus control (below 1 favors the checkpoint arm). For response endpoints put odds ratios in the
   hr_* columns and say so in notes. Leave a cell empty when not reported; never impute.
2. `screening.md`: the exact PubMed query run, the number of hits, how titles/abstracts were screened, every seed-list
   trial with its decision and reason, and any trials added from reference lists.
3. `summary.md`: ten lines at most: numbers screened, assessed, eligible (survival / response), and any doubts.

Work carefully rather than fast; a missed eligible trial or a mis-copied confidence interval matters more than time.
