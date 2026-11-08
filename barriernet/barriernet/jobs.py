"""Running a benchmark in its own process, for the web app: the page stays responsive, progress goes to a small text
file that the page reads every second, and the result comes back as a pickle.

The worker is a fresh ``python -m barriernet.jobs`` command (not multiprocessing's spawn, which re-imports the server's
main module): it is single-threaded when it forks its own workers, and it runs in its own process group, so cancelling
stops its workers too."""
from __future__ import annotations
import json
import os
import pathlib
import pickle
import shutil
import signal
import subprocess
import sys
import tempfile
from typing import Optional, Tuple


def _work(spec: dict, folder: pathlib.Path) -> None:
    def prog(stage, i, n):
        tmp = folder / "progress.tmp"; tmp.write_text(f"{stage}\t{i}\t{n}"); os.replace(tmp, folder / "progress.txt")
    try:
        from .benchmark import benchmark
        from .library import combine, load_library, read_signature_table
        prog("reading the data", 0, 1)
        lib = load_library()
        if spec.get("keep"):
            lib = lib[lib.signature.isin(spec["keep"])]
        if spec.get("signatures"):
            lib = combine(lib, read_signature_table(spec["signatures"]))
        if spec["source"] == "demo":
            from .demo import demo_trials
            trials = demo_trials()
        else:
            from .bundle import load_trials_zip
            trials, _ = load_trials_zip(spec["zip"])
        B = benchmark(trials, library=lib, n_random=spec["n_random"], n_boot=spec["n_boot"], seed=spec["seed"],
                      n_jobs=spec["n_jobs"], progress=prog, **({"genes": spec["genes"]} if spec.get("genes") else {}))
        tmp = folder / "result.tmp"
        with open(tmp, "wb") as f:
            pickle.dump(B, f)
        os.replace(tmp, folder / "result.pkl")
    except Exception as e:                                     # reported to the page
        (folder / "error.txt").write_text(f"{type(e).__name__}: {e}")


class Job:
    """A benchmark running as a separate Python process."""

    def __init__(self, spec: dict):
        self.spec = dict(spec)
        self.folder = pathlib.Path(tempfile.mkdtemp(prefix="barriernet_job_"))
        for key in ("zip", "signatures"):                     # uploads may be removed by the server once the request ends
            if self.spec.get(key):
                dst = self.folder / f"{key}{pathlib.Path(self.spec[key]).suffix}"; shutil.copy(self.spec[key], dst); self.spec[key] = str(dst)
        (self.folder / "spec.json").write_text(json.dumps(self.spec))
        root = str(pathlib.Path(__file__).resolve().parent.parent)          # the folder holding the barriernet package
        extra = [os.environ["PYTHONPATH"]] if os.environ.get("PYTHONPATH") else []
        env = dict(os.environ, PYTHONPATH=os.pathsep.join([root] + extra))
        self.log = open(self.folder / "worker.log", "wb")
        self.proc = subprocess.Popen([sys.executable, "-m", "barriernet.jobs", str(self.folder)], env=env, cwd=str(self.folder),
                                     stdout=self.log, stderr=subprocess.STDOUT, start_new_session=True)

    def progress(self) -> Optional[Tuple[str, int, int]]:
        try:
            s, i, n = (self.folder / "progress.txt").read_text().split("\t"); return s, int(i), int(n)
        except (OSError, ValueError):
            return None

    def done(self) -> bool:
        return self.proc.poll() is not None

    def result(self):
        self.proc.wait(timeout=5)
        if (self.folder / "error.txt").exists():
            raise RuntimeError((self.folder / "error.txt").read_text())
        if not (self.folder / "result.pkl").exists():
            log = self.folder / "worker.log"
            tail = log.read_text(errors="replace")[-400:] if log.exists() else ""
            raise RuntimeError("the benchmark stopped without a result. " + tail)
        with open(self.folder / "result.pkl", "rb") as f:
            return pickle.load(f)

    def cancel(self) -> None:
        if self.proc.poll() is None:
            try:
                os.killpg(self.proc.pid, signal.SIGTERM)
            except (OSError, AttributeError):
                self.proc.terminate()
            try:
                self.proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(self.proc.pid, signal.SIGKILL)
                except (OSError, AttributeError):
                    self.proc.kill()

    def cleanup(self) -> None:
        try:
            self.log.close()
        except Exception:
            pass
        shutil.rmtree(self.folder, ignore_errors=True)


if __name__ == "__main__":
    _folder = pathlib.Path(sys.argv[1])
    _work(json.loads((_folder / "spec.json").read_text()), _folder)
