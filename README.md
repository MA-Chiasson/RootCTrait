# RootCTrait

**Extraction of 3D root system architecture (RSA) traits from segmented CT
volumes, for population-scale phenotyping and quantitative genetics (GWAS).**

RootCTrait takes already segmented binary volumes of root systems (for example
from X-ray CT of plants grown in pots), reconstructs each skeleton, decomposes it
into ordered roots (primary root and laterals), detects the collar and removes the
hypocotyl, decontaminates surface artifacts, and writes an Excel table with 38
architectural traits per sample. Processing settings are shared across all
batches so that traits are directly comparable across a whole population.

The pipeline is **multi-batch**: one run processes one or several batches (groups
of samples), each in its own folder, and writes one Excel table per batch.

---

## Contents

- [Two ways to run RootCTrait](#two-ways-to-run-rootctrait)
- [Processing overview](#processing-overview)
- [Expected input](#expected-input)
- [Installation](#installation)
- [Quick start](#quick-start)
- [The desktop app (RootCTraitV3)](#the-desktop-app-rootctraitv3)
- [The command line](#the-command-line)
- [Configuration reference](#configuration-reference)
- [Worked examples](#worked-examples)
- [Output files](#output-files)
- [Performance and parallel processing](#performance-and-parallel-processing)
- [Reproducibility](#reproducibility)
- [Troubleshooting](#troubleshooting)
- [Notes for downstream analysis (GWAS)](#notes-for-downstream-analysis-gwas)
- [Repository structure](#repository-structure)
- [Known limitations](#known-limitations)
- [Credits, citation, license](#credits-and-acknowledgments)

---

## Two ways to run RootCTrait

RootCTrait has one engine (`run_pipeline.py` and the `rootctrait/` package) and
two front ends that share it:

1. **The desktop app, RootCTraitV3** (`rootctrait_app.py`): a graphical interface
   where you import your batch folders, set parameters, and watch a live console
   with a progress bar. Best for day-to-day work and visual review.
2. **The command line** (`python run_pipeline.py`): driven entirely by a
   `params.txt` file. Best for reproducible, scriptable runs on a server or an HPC
   node, and for the exact settings you report in a paper.

Both call the same processing code, so the traits they produce are identical for
identical settings. Pick whichever fits the moment; you can even start in the app
and finish on a server.

---

## Processing overview

For each sample:

1. **Load** the volume (multi-format, see below), binarize (`V > 0.5 * max(V)`)
   and crop to the bounding box. An empty mask (no voxel above the threshold) is
   reported as an error and skipped, so one bad file never stops a batch.
2. **Skeletonize** and prune short spurious branches.
3. **Detect the collar** at the thickest skeleton point within the uppermost 20%
   of the system in depth.
4. **Decompose** the skeleton into ordered segments (order 1 = primary root,
   order 2 = laterals, etc.).
5. **Decontaminate**: remove parallel sheets and floating fragments.
6. **Raise the collar and exclude the hypocotyl**: the collar climbs along the
   thick base column and stops at the hypocotyl; the hypocotyl (vertical stem
   column above the collar plus the branches hanging high on it) is excluded,
   while basal roots at collar level are kept.
7. **Orphan cleanup**: after the hypocotyl is removed, one connectivity pass drops
   any detached fragment. It is applied to every sample; its memory use grows
   linearly with the size of the system.
8. **Extract all traits** on the cleaned root skeleton, using the raised collar as
   the reference point; the primary root is traced continuously from the raised
   collar following the real skeleton path, and its whole length is measured. An
   upward return at its tip is not corrected: its height is reported in the
   `pivot_return` column and, above 3 mm, drawn in the review figure so that the
   sample can be checked. Root volume and surface
   are measured on the part of the mask attached to the cleaned roots.

A resume mechanism (`checkpoint_traits.jsonl`, per batch) lets you interrupt and
restart without recomputing everything.

---

## Expected input

RootCTrait expects an **already segmented binary foreground mask**, not a raw
grayscale CT volume. A multi-label segmentation must be reduced to the root label
beforehand.

Axis convention: `col0 = Y (depth)`, `col1 = X`, `col2 = Z`. The default voxel
size is `0.39, 0.39, 0.2` mm (configurable).

### Supported formats

The file extension determines the reader used. For `.mat` files holding several
variables, the largest 3D array is used.

| Extension         | Format                                | Library      |
|-------------------|---------------------------------------|--------------|
| `.mat`            | MATLAB v7.3 (HDF5) or older (5/6/7)   | h5py / scipy |
| `.tif`, `.tiff`   | 3D TIFF stack                         | tifffile     |
| `.npy`, `.npz`    | NumPy array                           | numpy        |
| `.nii`, `.nii.gz` | NIfTI                                 | nibabel      |

`tifffile` and `nibabel` are imported only if the corresponding format is used.

### Axis orientation

If a loaded volume is not in the order `col0 = depth`, set the optional axis order
(for example `2,1,0`) on the batch line in `params.txt`, then check on the produced
HTML figure that depth points downward.

---

## Installation

```bash
git clone https://github.com/MA-Chiasson/RootCTrait.git
cd RootCTrait
pip install -r requirements.txt        # simplest: just the dependencies
```

Python 3.10 or newer is recommended. If you also want to use RootCTrait as a
library from your own scripts (`import rootctrait`), install the package instead:

```bash
pip install -e .        # dependencies + the importable rootctrait package
```

Dependencies: numpy, scipy, scikit-image, h5py, pandas, openpyxl, plotly
(see `requirements.txt`). The desktop app additionally uses Tkinter, which ships
with standard Python, and installs `tkinterdnd2` on first launch for drag-and-drop.

---

## Quick start

### With the desktop app

```bash
python rootctrait_app.py
```

Then: **Import Folders** to add your batch folders, **Parameters** to choose the
output folder, and **Analyse** to run. Details in
[The desktop app](#the-desktop-app-rootctraitv3).

### With the command line

```bash
python run_pipeline.py            # reads params.txt in the same folder
```

Point it at a different config with the `PARAMS` environment variable (here the
bundled synthetic example, which runs out of the box):

```bash
PARAMS=example/params_example.txt python run_pipeline.py
```

---

## The desktop app (RootCTraitV3)

Launch it with `python rootctrait_app.py`. The window has two tabs.

### Analysis tab

The workflow reads left to right:

- **New Session**: clears the queue, the console and the saved workspace, and
  forgets the last parent folder, so the next Import Folders opens empty.
- **Import Folders**: opens a picker. Click **Browse Parent Folder** to choose the
  folder that contains your batches (for example `data`), then tick the batch
  subfolders you want in the queue. Each ticked folder becomes one batch, named
  after the folder.
- **Parameters**: opens the settings window (see below). You must choose an
  **output folder** here before you can analyse; until you do, the status reads
  `Set output folder in Parameters`.
- **Analyse / Stop**: start or stop processing. Progress, a percentage and an ETA
  appear on the progress bar; every sample prints a line in the live console.
- **Save Console Log**: writes the console text to a file.

The queue table shows each batch, its folder path and its status. The session
(queue, parameters, last parent folder) is saved automatically and restored the
next time you open the app.

### Parameters window

Fields marked `*` are required. The one you must set is the **output folder
(results root)**; the rest have sensible defaults.

| Field                        | Meaning                                             |
|------------------------------|-----------------------------------------------------|
| Output folder (results root) | where `results/<batch>/` is written (required)      |
| Voxel size (depth,x,z) mm    | physical voxel size, order Y,X,Z                    |
| File name pattern            | how sample files are named, with `{name}` as the id |
| Prune length (voxels)        | skeleton pruning length                             |
| Min segment length (mm)      | shortest segment kept                               |
| Sheet: parallel neighbours   | decontamination `BC_MIN`                            |
| Sheet: max linearity         | decontamination `LIN_MAX`                           |
| Sheet: max length (mm)       | decontamination `LEN_MAX`                           |
| Timeout per sample (s)       | per-sample time limit                               |
| Drop orphan fragments        | remove detached fragments (on/off)                  |
| Save 3D figures              | write the interactive HTML figures (on/off)         |

You can save the current settings as a named **parameter set** and reload it
later (the output folder is deliberately not stored in a preset, since it changes
from run to run). `PARALLEL` has no dedicated field; the app inherits it from
`params.txt` at startup, so set it there if you need it.

### Tools tab

Utilities that act on an existing `results` folder. The **Table format** selector
(both / xlsx / csv) applies to the tools that write tables. Hover any button for a
short description.

- **Merge batches**: merges the trait tables of all batches into one combined file.
- **Extract checkpoints**: rebuilds trait tables straight from the
  `checkpoint_traits.jsonl` checkpoints, without rerunning the analysis.
- **Figure report**: builds an HTML report of the interactive 3D Plotly figures,
  sorted by batch and in numerical order, for visual validation.
- **Delete checkpoints**: deletes the `checkpoint_traits.jsonl` checkpoints of a
  results folder to start over; progress is reset to 0 and the analysis restarts from the beginning.

---

## The command line

`python run_pipeline.py` reads all its settings from `params.txt` in the same
folder (or the file named by the `PARAMS` environment variable). The file has two
kinds of lines: `KEY=value` settings, and one `BATCH` (or `BLOC`) line per batch.

```bash
python run_pipeline.py                              # uses params.txt
PARAMS=example/params_example.txt python run_pipeline.py   # uses a named config
```

Nothing else is needed; the folders in `results/` are created automatically.

---

## Configuration reference

All settings live in `params.txt`. Keys are case-insensitive. Lines starting with
`;` are comments, and ` ; text` after a setting is ignored.

| Key             | Default        | Role                                                        |
|-----------------|----------------|------------------------------------------------------------|
| `BATCHES`       | `ALL`          | `ALL`, or a comma-separated list of batch names to process |
| `DATA_ROOT`     | `data`         | root folder holding `data/<batch>/`                        |
| `RESULTS_ROOT`  | `results`      | root folder for `results/<batch>/`                         |
| `VOXEL`         | `0.39,0.39,0.2`| voxel size in mm, order Y (depth), X, Z                     |
| `PRUNE_VOX`     | `5`            | skeleton pruning length (voxels)                           |
| `MIN_SEG_LEN_MM`| `2.0`          | shortest segment kept (mm)                                 |
| `BC_MIN`        | `3`            | decontamination: min parallel neighbours for a sheet       |
| `LIN_MAX`       | `0.7`          | decontamination: max linearity for a sheet                 |
| `LEN_MAX`       | `15`           | decontamination: max length for a sheet (mm)               |
| `DROP_ORPHANS`  | `1`            | drop floating fragments (1/0)                              |
| `SAVE_FIGURES`  | `1`            | write the interactive HTML figures (1/0)                   |
| `TIMEOUT`       | `1800`         | per-sample time limit (seconds)                            |
| `PARALLEL`      | `1`            | number of samples processed concurrently (1 = serial)      |

The processing settings above the batch lines are **shared by every batch**, which
is what keeps traits comparable across the whole population for GWAS.

### Defining batches

One line per batch, with the `BATCH` keyword (or its French alias `BLOC`):

```
BATCH <name> | <file_pattern> | <axis_order optional>
```

- `<name>` must match the subfolder in `data/` exactly.
- `{name}` in the pattern is the sample id; the extension sets the format.
- `<axis_order>` is optional (for example `2,1,0`); leave it empty for no
  permutation.

Example (`BATCH` and `BLOC` are interchangeable):

```
BATCH block1_t1 | block1_t1_{name}.mat |
BATCH block2_t1 | block2_t1_{name}.mat |
BATCH block3_t2 | block3_t2_{name}.tif | 2,1,0
```

Here `block1_t1_{name}.mat` matches `block1_t1_S001.mat`, `block1_t1_S002.mat`,
and so on, and the sample IDs are `S001`, `S002`, ...

Choose which of the defined batches to run with `BATCHES`:

```
; every batch defined below
BATCHES=ALL
; only these two
BATCHES=block1_t1,block3_t2
```

---

## Worked examples

### Example 1: one batch, from the command line

Layout:

```
data/
  mybatch/
    scan_S1.mat
    scan_S2.mat
    ...
```

`params.txt` (only the lines that matter here):

```
DATA_ROOT=data
RESULTS_ROOT=results
BATCHES=ALL
VOXEL=0.39,0.39,0.2
BATCH mybatch | scan_{name}.mat |
```

```bash
python run_pipeline.py
```

Output: `results/mybatch/traits_mybatch.xlsx` and `results/mybatch/figures/`.

### Example 2: several batches at once

One subfolder per batch, one batch line each; a single run processes them all with
identical settings, so the traits stay comparable across batches:

```
data/
  block1_t1/  block1_t1_S001.mat block1_t1_S002.mat ...
  block2_t1/  block2_t1_S001.mat block2_t1_S002.mat ...
  block3_t1/  block3_t1_S001.mat block3_t1_S002.mat ...
```

```
BATCHES=ALL
BATCH block1_t1 | block1_t1_{name}.mat |
BATCH block2_t1 | block2_t1_{name}.mat |
BATCH block3_t1 | block3_t1_{name}.mat |
```

```bash
python run_pipeline.py
```

Output: `results/block1_t1/`, `results/block2_t1/`, `results/block3_t1/`, each with
its own trait table, figures, and a `params_used.json`. Set
`BATCHES=block1_t1,block3_t1` to process only some of them.

### Example 3: run several batches faster (parallel)

Process four samples at a time. Every sample still runs in its own process, so the
per-sample timeout still applies:

```
PARALLEL=4
BATCHES=ALL
BATCH block1_t2 | block1_t2_{name}.mat |
BATCH block2_t2 | block2_t2_{name}.mat |
BATCH block3_t2 | block3_t2_{name}.mat |
```

```bash
python run_pipeline.py
```

Validate `PARALLEL` on one batch first, then scale up. See
[Performance](#performance-and-parallel-processing).

### Example 4: resume after an interruption

If a run is stopped (Ctrl-C, a crash, or **Stop** in the app), just launch it
again. Samples already in `checkpoint_traits.jsonl` are skipped and only the rest
are processed:

```bash
python run_pipeline.py        # picks up where it left off
```

To force a full recompute of a batch, delete its `checkpoint_traits.jsonl` (or use
**Delete checkpoints** in the app), then run again.

### Example 5: data and results anywhere on disk

Absolute paths are allowed, so the code folder can stay clean:

```
DATA_ROOT=/mnt/scans/soybean_masks
RESULTS_ROOT=/mnt/analysis/rootctrait_out
```

### Example 6: the Tools, from the command line

The same utilities exposed in the app's Tools tab can run standalone:

```bash
python -m tools.figure_report               # writes results/figure_report.html
python -m tools.merge_batches --format both # merges all batch tables into one file
python -m tools.extract_checkpoint          # rebuilds tables from the checkpoints
```

### Example 7: sensitivity of the traits to the thresholds

`tools/sensitivity.py` moves each threshold below and above its default (one at a
time), reruns the pipeline on a random subset of samples, and reports for every
trait the rank correlation with the default and the median relative change:

```bash
python -m tools.sensitivity --batch block1_t1 --batch block2_t2 --n 30 --workers 4
```

Output: `results/sensitivity/sensitivity_traits.csv` and
`results/sensitivity/sensitivity_summary.csv`.

---

## Output files

For each batch, in `results/<batch>/`:

- **`traits_<batch>.xlsx`**: trait table, one row per sample (lengths in cm,
  diameters in mm, volumes in cm3, angles in degrees). Columns `n_raw`
  (segments before cleaning), `n_removed` and `%removed` report the
  decontamination, and `pivot_return` (mm) the upward return of the pivot after
  its deepest point (0 when the deepest point is the tip). Tables and checkpoints
  written by versions before 2.1.0, which
  used French column names (`n_brut`, `n_retire`, `%retire`, `NRL_court_<5`,
  `NRL_moyen_5_15`), are still read and renamed by the tools.
- **`figures/<sample>.html`**: interactive 3D view. Pivot (black, from the raised
  collar, as used for the traits), upward return of the pivot tip above 3 mm, kept
  in the traits but flagged for checking (brown, dotted), kept laterals (blue), removed pollution (red), hypocotyl (orange,
  excluded), detached orphans (grey), original collar (green), raised collar
  (purple diamond).
- **`checkpoint_traits.jsonl`**: resume state (delete to recompute).
- **`params_used.json`**: the exact settings used for this batch (voxel, pruning,
  decontamination thresholds, timeout, parallel, orphan threshold, data folder,
  timestamp). Keep it with your results for reproducibility and the methods
  section of a paper.
- **`failures.jsonl`**: one line per sample that timed out or errored, with the
  name, the failure type, the message and a timestamp. Empty or absent when
  everything succeeds; a quick way to see what needs attention.

See [`docs/traits.md`](docs/traits.md) for the full trait list with definitions
and units.

---

## Performance and parallel processing

By default `PARALLEL=1`: samples are processed one at a time. This is the safe,
proven mode and the recommended default for a final production run.

Set `PARALLEL=N` to process N samples concurrently. The scheduler keeps at most N
processes alive, and because each sample still runs in its own process, the
per-sample timeout-kill is preserved and the checkpoint is written by the single
parent process (no locking, no corruption). A good starting point is the number of
physical cores minus one.

Guidance:

- **Validate first.** Run one batch with your chosen `PARALLEL` and confirm the
  trait table matches a serial run before relying on it at scale.
- Serial and parallel produce the **same set of results**; only the wall-clock
  time differs.
- Very large single systems are memory-heavy; if you see memory pressure, lower
  `PARALLEL` or raise per-sample limits.

---

## Reproducibility

Every batch writes a `params_used.json` capturing the effective settings for that
run. Because all processing parameters travel with each sample to its worker
process, the settings are honored identically whether Python uses the `fork`
(Linux) or `spawn` (Windows, macOS) start method: what you set in the app or in
`params.txt` is exactly what the computation uses.

For a citable, frozen configuration, keep the `params.txt` (or the relevant
`params_*.txt`) alongside your results, together with the `params_used.json` files.

---

## Troubleshooting

- **"No batch to process."** The `BATCHES` list is empty, or no `BATCH`/`BLOC`
  line matched. Check that each batch name matches a subfolder in `DATA_ROOT`, and
  that the batch lines use `BATCH` or `BLOC` (both are accepted).
- **"no file matching ... "** The `{name}` pattern or the extension does not match
  the files in the batch folder. Check the pattern against a real filename.
- **`empty mask` error for a sample.** The file has no voxel above the binarization
  threshold (empty or non-binary input). It is skipped and logged in
  `failures.jsonl`; the rest of the batch continues.
- **Timeouts.** Increase `TIMEOUT`, or lower `PARALLEL` if the machine is
  oversubscribed. Timed-out samples are retried on the next run and recorded in
  `failures.jsonl`.
- **Figures not written.** `SAVE_FIGURES=0`; set it to `1` (or tick "Save 3D
  figures" in the app). Plotly must be installed.
- **Depth looks inverted in a figure.** Set the batch `axis_order` (for example
  `2,1,0`) and check again.

---

## Notes for downstream analysis (GWAS)

- Trait quality depends on segmentation quality. Traits that **aggregate over all
  segments** (total length, counts, branching order, density) get inflated by
  residual surface pollution; traits that measure a **specific geometry** (primary
  root length, branching angles) are more robust. Volume traits (root volume, convex
  hull, surface, compactness) are the most sensitive.
- The size and count traits are **highly redundant** (they mostly measure one
  "system size" axis). Prefer a small non-redundant set over the full list.
- Treat `%removed`, `MaxO` and collar-related quantities as covariates or quality
  indicators rather than biological traits. Including the batch and `%removed` as
  covariates removes most of the scan-quality confound.
- Keep `params_used.json` with your trait tables so the exact processing settings
  are traceable per batch.

---

## Repository structure

```
.
├── rootctrait_app.py           Desktop app (RootCTraitV3): import, parameters, run, tools
├── pipeline_api.py             Callable API layer used by the app (and any script)
├── run_pipeline.py             Engine: I/O, collar, hypocotyl, traits, Excel, batch loop
├── rootctrait/                 Pipeline package (pip-installable, `import rootctrait`)
│   ├── __init__.py
│   ├── io_volume.py                Multi-format loading of 3D volumes
│   ├── graph_extraction.py         Skeleton graph, branch points, pruning
│   ├── root_decomposition.py       Decomposition into ordered roots
│   ├── decontamination.py          Parallel sheets + orphan fragments
│   ├── detection_hypocotyle.py     Bounded collar + hypocotyl detection
│   ├── root_traits_full.py         Full trait set
│   └── legacy.py                   Former column names (reads files from versions < 2.1.0)
├── tools/                      Post-processing utilities (also in the app's Tools tab)
│   ├── __init__.py
│   ├── merge_batches.py            Merge all batch tables into one file
│   ├── extract_checkpoint.py       Rebuild tables from checkpoints
│   ├── figure_report.py            HTML index to review 3D figures (QC)
│   └── sensitivity.py              One-at-a-time sensitivity analysis of the thresholds
├── docs/traits.md              Trait reference
├── docs/limitations.md         Known limitations
├── params.txt                  All settings + batch definitions
├── example/params_example.txt  Ready-to-run config for the bundled example
├── example/roots/sample_S1.npy Synthetic example dataset (versioned)
├── tests/test_invariants.py    Invariant tests (pytest)
├── tests/test_components.py    Unit tests: connectivity, pivot hook, root mask
├── validation/validate_phantoms.py  Accuracy check on known-geometry phantoms
├── validation/phantom_results.csv   Results of that check
├── pyproject.toml              Package metadata (pip install -e .)
├── requirements.txt
├── CITATION.cff                Citation metadata
├── data/                       Input volumes (not versioned)
└── results/                    Output: Excel, figures, checkpoint, params_used.json, failures.jsonl
```

Run `run_pipeline.py` or `rootctrait_app.py` from the root folder (or
`pip install -e .` first) so that `rootctrait` and `tools` are importable.

---

## Known limitations

RootCTrait is a research tool with clearly stated limits: it depends entirely on
the upstream segmentation, oblique or poorly segmented hypocotyls can be missed,
thresholds are calibrated on soybean CT at ~0.39 x 0.39 x 0.2 mm, and there is no
ground truth for 3D roots (validation is by consistency, reproducibility and
visual inspection, not absolute accuracy). See [`docs/limitations.md`](docs/limitations.md)
for the full discussion.

---

## Credits and acknowledgments

This pipeline builds on the work of **Mana Eskandari** on root system architecture
phenotyping, whose skeleton-reconstruction foundations (`graph_extraction.py`) are
reused here.

The related project **RootWeave** (Xuehai Zhou et al., *Computers and Electronics
in Agriculture*, 2025, https://github.com/xuehai-zhou/RootWeave) is a useful
reference and a future integration path.

Developed as part of doctoral work at Université Laval (CT analysis of soybean
root systems).

## Citation

If you use RootCTrait, please cite the accompanying article (in preparation) and
the software itself. A machine-readable citation is in [`CITATION.cff`](CITATION.cff).

Each release is archived on Zenodo with a permanent DOI: version 2.0.0 is
[10.5281/zenodo.23019246](https://doi.org/10.5281/zenodo.23019246) and version
1.0.0 is [10.5281/zenodo.22212849](https://doi.org/10.5281/zenodo.22212849). Cite the
version you used, for example: "Root traits were extracted with RootCTrait v2.0.0
(doi:10.5281/zenodo.23019246)."

## License

Released under the MIT License. See [`LICENSE`](LICENSE). Copyright (c) 2026
Marc-Antoine Chiasson.
