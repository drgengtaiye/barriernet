"""Reading a zip of trials or cohorts (a YAML file plus the tables it names) for the web app."""
from __future__ import annotations
import pathlib
import tempfile
import zipfile
import os
from typing import List, Tuple

MAX_MEMBERS = int(os.environ.get("BARRIERNET_MAX_ZIP_FILES", 2000))
MAX_UNZIPPED_MB = int(os.environ.get("BARRIERNET_MAX_UNZIPPED_MB", 2000))


def extract(zip_path, dest=None) -> pathlib.Path:
    """Extract into a new folder, refusing members whose path would leave it."""
    dest = pathlib.Path(dest or tempfile.mkdtemp(prefix="barriernet_upload_")).resolve()
    with zipfile.ZipFile(zip_path) as z:
        members = z.infolist()
        if len(members) > MAX_MEMBERS or sum(m.file_size for m in members) > MAX_UNZIPPED_MB * 1024 * 1024:
            raise ValueError(f"the zip holds more than {MAX_MEMBERS} files or unpacks to more than {MAX_UNZIPPED_MB} MB")
        for m in members:
            target = (dest / m.filename).resolve()
            if dest not in target.parents and target != dest:
                raise ValueError(f"unsafe path in zip: {m.filename}")
        z.extractall(dest)
    return dest


def _config(folder: pathlib.Path, key: str) -> Tuple[dict, pathlib.Path]:
    import yaml
    for y in sorted(folder.rglob("*.y*ml")):
        cfg = yaml.safe_load(open(y, encoding="utf-8"))
        if isinstance(cfg, dict) and key in cfg:
            base = y.parent
            for item in cfg[key]:
                for f in ("expr", "clinical"):
                    p = (base / str(item[f])).resolve()
                    if folder not in p.parents:
                        raise ValueError(f"{item.get('name')}: {f} points outside the uploaded folder")
                    if not p.exists():
                        raise ValueError(f"{item.get('name')}: {item[f]} not found in the zip")
                    item[f] = str(p)
            return cfg, y
    raise ValueError(f"no YAML file with a '{key}:' list in the zip")


def _load(zip_path, key, loader):
    import shutil
    folder = extract(zip_path)
    try:
        cfg, _ = _config(folder, key)
        return loader(cfg), cfg
    finally:                                                   # the tables are in memory now
        shutil.rmtree(folder, ignore_errors=True)


def load_trials_zip(zip_path) -> Tuple[list, dict]:
    from .cli import _load_trials
    return _load(zip_path, "trials", _load_trials)


def load_cohorts_zip(zip_path) -> Tuple[List, dict]:
    from .cli import _load_cohorts
    return _load(zip_path, "cohorts", _load_cohorts)
