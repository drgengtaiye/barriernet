# barriernet web app

A Shiny for Python app (`app.py`) for the analyses of barriernet 0.4: the benchmark of published signatures in
randomized trials (synthetic demonstration or your own trials), single-arm cohorts, pooling of published subgroup
results with an editable table, the trial-design calculator and the signature library. No R packages are used.

It runs in the browser at https://drgengtaiye.github.io/barriernet/ (built by `.github/workflows/pages.yml` with Shinylive, Python compiled
to WebAssembly). There nothing leaves the visitor's computer; a browser has no worker processes, so the benchmark shows
a demonstration computed in advance at once and runs uploaded data in the page on one core, with at most 200 random sets
and 100 resamples.

Run it on your computer from this folder:

```bash
python -m pip install -r requirements.txt
shiny run app.py                 # or, with the package installed: barriernet app
```

For Shiny Server (Open Source 1.5.21 or later hosts Shiny for Python apps), make a self-contained folder and copy it to
the server's app directory:

```bash
./make_bundle.sh /path/to/bundle
```

The server needs a Python with the packages in `requirements.txt`, set for this location with Shiny Server's `python`
setting. Everything the app reads is inside the bundle (the `barriernet` package and its `barriernet/data/` tables,
read by relative paths). It writes only temporary files under the system temporary directory and removes them when a
run or a browser session ends. The benchmark runs as a separate `python -m barriernet.jobs` process.

Limits for a shared server, each changeable by an environment variable: `BARRIERNET_MAX_WORKERS` (worker processes per
benchmark; default min(4, CPUs − 1)), `BARRIERNET_MAX_JOBS` (benchmarks at once; 1), `BARRIERNET_MAX_RANDOM` (1000),
`BARRIERNET_MAX_BOOT` (300), `BARRIERNET_MAX_UPLOAD_MB` (200), `BARRIERNET_MAX_UNZIPPED_MB` (2000),
`BARRIERNET_MAX_ZIP_FILES` (2000).

Memory: loading the app takes about 235 MB; a benchmark of the synthetic demonstration with four workers peaks near
850 MB summed over its processes. Larger uploaded trials need more.
