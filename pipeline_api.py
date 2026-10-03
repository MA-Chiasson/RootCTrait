"""pipeline_api.py

Callable API layer over run_pipeline, for the desktop GUI (and any other caller).

It does not modify the command-line script: the CLI (python run_pipeline.py) keeps
working unchanged. This layer lets a caller:
  - read the default parameters (default_params),
  - run an analysis on imported folders, with parameters passed in,
  - receive progress messages live through a callback (so the GUI can show them
    in its terminal pane), instead of only printing to the console.

Each imported folder is treated as one batch; the folder name becomes the batch
name, verbatim (no naming convention assumed).
"""
import os
import contextlib
import run_pipeline as rp


def default_params():
    """Current default parameters (from params.txt if present, else built-in)."""
    return {
        "VOXEL": ",".join(str(x) for x in rp.VOXEL_SIZE),
        "PRUNE_VOX": rp.PRUNE_VOX,
        "MIN_SEG_LEN_MM": rp.MIN_SEG_LEN_MM,
        "BC_MIN": rp.BC_MIN,
        "LIN_MAX": rp.LIN_MAX,
        "LEN_MAX": rp.LEN_MAX,
        "DROP_ORPHANS": rp.DROP_ORPHANS,
        "DENS_MAX": rp.DENS_MAX,
        "DENS_LEN_MAX": rp.DENS_LEN_MAX,
        "RESCUE_MIN_MM": rp.RESCUE_MIN_MM,
        "SAVE_FIGURES": rp.SAVE_FIGURES,
        "TIMEOUT": rp.TIMEOUT,
        "PATTERN": "{name}.mat",
        "DATA_ROOT": rp.DATA_ROOT,
        "RESULTS_ROOT": rp.RESULTS_ROOT,
        "PARALLEL": rp.PARALLEL,
    }


class _LineWriter:
    """A file-like object that forwards complete lines to a callback."""
    def __init__(self, cb):
        self.cb = cb
        self.buf = ""

    def write(self, s):
        self.buf += s
        while "\n" in self.buf:
            line, self.buf = self.buf.split("\n", 1)
            self.cb(line)

    def flush(self):
        if self.buf:
            self.cb(self.buf)
            self.buf = ""


def _apply(params):
    """Push the parameter dict into run_pipeline's module globals."""
    rp.VOXEL_SIZE = tuple(float(x) for x in str(params.get("VOXEL", "0.39,0.39,0.2")).split(","))
    rp.PRUNE_VOX = int(params.get("PRUNE_VOX", 5))
    rp.MIN_SEG_LEN_MM = float(params.get("MIN_SEG_LEN_MM", 2.0))
    rp.BC_MIN = int(params.get("BC_MIN", 3))
    rp.LIN_MAX = float(params.get("LIN_MAX", 0.7))
    rp.LEN_MAX = float(params.get("LEN_MAX", 15))
    rp.DROP_ORPHANS = bool(params.get("DROP_ORPHANS", True))
    rp.DENS_MAX = float(params.get("DENS_MAX", 35))
    rp.DENS_LEN_MAX = float(params.get("DENS_LEN_MAX", 6))
    rp.RESCUE_MIN_MM = float(params.get("RESCUE_MIN_MM", 5))
    rp.SAVE_FIGURES = bool(params.get("SAVE_FIGURES", True))
    rp.TIMEOUT = int(params.get("TIMEOUT", 1800))
    rp.PARALLEL = max(1, int(params.get("PARALLEL", 1)))
    rp.DATA_ROOT = params.get("DATA_ROOT", "data") or "data"
    rp.RESULTS_ROOT = params.get("RESULTS_ROOT", "results")


def run(folders, params=None, progress=None, should_stop=None):
    """Run the pipeline on the given folders.

    folders  : list of folder paths; each folder = one batch, name = folder name
    params   : parameter dict (see default_params); None uses defaults
    progress : callback taking one string (a progress line); None prints to console
    should_stop : callback returning True if stop requested; None ignores stops
    """
    params = params or default_params()
    progress = progress or (lambda s: print(s))
    should_stop = should_stop or (lambda: False)

    _apply(params)
    pattern = params.get("PATTERN", "{name}.mat")
    writer = _LineWriter(progress)
    with contextlib.redirect_stdout(writer):
        if not folders:
            print("No folder to process.")
        for folder in folders:
            if should_stop():
                print("Analysis cancelled by user.")
                break
            folder = os.path.abspath(folder)
            rp.DATA_ROOT = os.path.dirname(folder)          # parent -> DATA_ROOT/name = folder
            batch = {"name": os.path.basename(folder),
                     "pattern": pattern, "axis_order": None}
            rp.process_batch(batch, should_stop=should_stop)
        print("Done.")
    writer.flush()
