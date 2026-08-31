# RootCTrait

**Extraction of 3D root system architecture (RSA) traits from segmented CT
volumes, for population-scale phenotyping and quantitative genetics (GWAS).**

RootCTrait takes already segmented binary volumes of root systems (for example
from X-ray CT of plants grown in pots), reconstructs each skeleton, decomposes it
into ordered roots (primary root and laterals), detects the collar and removes the
hypocotyl, decontaminates surface artifacts, and writes an Excel table with about
40 architectural traits per sample. Processing settings are shared across all
batches so that traits are directly comparable across a whole population.

The pipeline is **multi-batch**: one run processes one or several batches (groups
of samples), each in its own folder, and writes one Excel table per batch.

---

## Processing overview

For each sample:

1. **Load** the volume (multi-format, see below), binarize (`V > 0.5 * max(V)`)
   and crop to the bounding box.
2. **Skeletonize** and prune short spurious branches.
3. **Detect the collar** at the thickest point of the upper region, then **raise
   it** along the thick base column, stopping at the hypocotyl.
4. **Decompose** the skeleton into ordered segments (order 1 = primary root,
   order 2 = laterals, etc.).
5. **Decontaminate**: remove parallel sheets and floating fragments.
6. **Identify and exclude the hypocotyl** (vertical stem column above the collar
   plus the branches hanging high on it); basal roots at collar level are kept.
7. **Orphan cleanup**: after the hypocotyl is removed, one connectivity pass drops
   any detached fragment.
8. **Extract all traits** on the cleaned root skeleton, using the raised collar as
   the reference point; the primary root is traced continuously from the raised
   collar following the real skeleton path.

A resume mechanism (`checkpoint_traits.jsonl`, per batch) lets you interrupt and
restart without recomputing everything.

## Expected input

RootCTrait expects an **already segmented binary foreground mask**, not a raw
grayscale CT volume. A multi-label segmentation must be reduced to the root label
beforehand.

Axis convention: `col0 = Y (depth)`, `col1 = X`, `col2 = Z`. The default voxel
size is `0.39, 0.39, 0.2` mm (configurable in `params.txt`).

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
on the batch line in `params.txt` (for example `2,1,0`), then check on the produced
HTML figure that depth points downward.

## Installation

```bash
git clone https://github.com/MA-Chiasson/RootCTrait.git
cd RootCTrait
pip install -r requirements.txt        # simplest: just the dependencies
```

Then run the pipeline from the repo folder with `python run_pipeline.py`.

Python 3.10 or newer recommended. If you also want to use RootCTrait as a library
from your own scripts (`import rootctrait`), install the package instead:

```bash
pip install -e .        # dependencies + the importable rootctrait package
```

## Try it on the bundled example

A small synthetic root system is included so you can run the whole pipeline in a
few seconds without any data of your own:

```bash
python run_pipeline.py            # uses params.txt
# or, explicitly, the bundled example:
PARAMS=params_example.txt python run_pipeline.py
```

This processes `example/roots/sample_S1.npy` and writes
`results/roots/traits_roots.xlsx` plus a 3D figure.

## Tests

Invariant tests run the pipeline on the bundled example and check the
internal-consistency rules (TRL >= LRP, IC <= 1, NRL = sum of length classes,
angles in [0, 90], no NaN in core traits):

```bash
pytest -q            # or: python tests/test_invariants.py
```

## Folder layout

```
data/<batch>/       the masks of this batch
results/<batch>/    output for this batch (created automatically):
                    traits_<batch>.xlsx, figures/, checkpoint_traits.jsonl
```

Run `bash makeDir.sh` to create `data/` and `results/`.

## Configuration

All settings are in `params.txt`, read by `run_pipeline.py` (same folder).

| Key                            | Role                                                       |
|--------------------------------|------------------------------------------------------------|
| `BATCHES`                      | `ALL`, or a comma-separated list of batch names            |
| `DATA_ROOT`, `RESULTS_ROOT`    | root folders for input and output                          |
| `VOXEL`                        | voxel size in mm, order Y,X,Z                              |
| `PRUNE_VOX`, `MIN_SEG_LEN_MM`  | skeleton pruning and minimum segment length               |
| `BC_MIN`, `LIN_MAX`, `LEN_MAX` | decontamination rule parameters                            |
| `DROP_ORPHANS`, `SAVE_FIGURES` | remove floating fragments; write HTML figures (1/0)        |
| `TIMEOUT`                      | per-sample time limit (seconds)                            |

Each batch is one line:

```
BATCH <name> | <file_pattern> | <axis_order optional>
```

`<name>` must match the subfolder in `data/`; `{name}` in the pattern is the
sample id; the extension sets the format. Example:

```
BATCH batch1 | sample_{name}.mat |
BATCH batch2 | scan{name}.tif    |
```

Processing settings above the batch lines are shared by every batch, which keeps
traits comparable across the whole population.

## Usage

Put your masks in `data/<batch>/`, edit `params.txt`, then run
`python run_pipeline.py`. Three worked examples:

### Example 1: the bundled synthetic example

Runs out of the box, no data of your own needed:

```bash
PARAMS=params_example.txt python run_pipeline.py
```

Output goes to `results/roots/traits_roots.xlsx` plus a 3D figure.

### Example 2: one batch of your own scans

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
BATCH mybatch | scan_{name}.mat |
```

```bash
python run_pipeline.py
```

Output: `results/mybatch/traits_mybatch.xlsx` and `results/mybatch/figures/`.

### Example 3: several batches at once

One subfolder per batch, one `BATCH` line each; a single run processes them all
with identical settings (so the traits stay comparable across batches):

```
data/
  march/   scan_S1.mat scan_S2.mat ...
  april/   scan_S1.mat scan_S2.mat ...
  june/    scan_S1.mat scan_S2.mat ...
```

```
BATCHES=ALL
BATCH march | scan_{name}.mat |
BATCH april | scan_{name}.mat |
BATCH june  | scan_{name}.mat |
```

```bash
python run_pipeline.py
```

Output: `results/march/`, `results/april/`, `results/june/`, each with its own
trait table and figures. Set `BATCHES=march,june` to process only some of them.

To review all the figures across batches and rate them (Good / Doubtful / Bad,
with a CSV export):

```bash
python generer_rapport_figures.py   # writes results/rapport_figures.html
```

### Keeping your folder tidy

Only `data/` and `results/` grow as you use the tool, and both are ignored by
git, so a clone stays clean. If you would rather keep the code folder untouched,
point `DATA_ROOT` and `RESULTS_ROOT` at folders anywhere on your machine
(absolute paths are allowed).

## Output

For each batch, in `results/<batch>/`:

- `traits_<batch>.xlsx`: trait table, one row per sample (lengths in cm, diameters
  in mm, volumes in cm3, angles in degrees).
- `figures/<sample>.html`: interactive 3D view. Pivot (black, from the raised
  collar), kept laterals (blue), removed pollution (red), hypocotyl (orange,
  excluded), detached orphans (grey), original collar (green) and raised collar
  (purple diamond).
- `checkpoint_traits.jsonl`: resume state (delete to recompute).

See [`docs/traits.md`](docs/traits.md) for the full trait list with definitions
and units.

### Visual quality control

`generer_rapport_figures.py` builds a single `rapport_figures.html` index to review
every 3D figure without loading them all at once, with Good / Doubtful / Bad
buttons and a free note per sample, and a CSV export of the judgments. Run it from
the root after processing:

```bash
python generer_rapport_figures.py
```

## Validation on synthetic phantoms

Because there is no ground truth for real 3D roots, accuracy is assessed on
synthetic root phantoms of known geometry (known total length, number of laterals,
and branching angles). The pipeline is run on each phantom and the measured traits
are compared to the true values:

```bash
python validation/validate_phantoms.py    # writes validation/phantom_results.csv
```

On the six bundled phantoms the mean absolute error is about **2% for primary root
length, 2% for total root length, 0% for the number of laterals (recovered exactly),
and 5% for the order-2 branching angle**. This measures the recovery accuracy of the
pipeline itself, independently of any upstream segmentation.

## Notes for downstream analysis (GWAS)

- Trait quality depends on segmentation quality. Traits that **aggregate over all
  segments** (total length, counts, branching order, density) get inflated by
  residual surface pollution; traits that measure a **specific geometry** (primary
  root length, branching angles) are more robust. Volume traits computed on the raw
  mask (root volume, convex hull, surface, compactness) are the most sensitive.
- The size and count traits are **highly redundant** (they mostly measure one
  "system size" axis). Prefer a small non-redundant set over the full list.
- Treat `%removed`, `MaxO` and collar-related quantities as covariates / quality
  indicators rather than biological traits. Including the batch and `%removed` as
  covariates removes most of the scan-quality confound.

## Repository structure

```
.
├── run_pipeline.py             Orchestration: I/O, collar, hypocotyl, traits, Excel
├── generer_rapport_figures.py  Builds an HTML index to review figures (QC)
├── rootctrait/                Pipeline package (pip-installable, `import rootctrait`)
│   ├── __init__.py
│   ├── io_volume.py                Multi-format loading of 3D volumes
│   ├── graph_extraction.py         Skeleton graph, branch points, pruning
│   ├── root_decomposition.py       Decomposition into ordered roots
│   ├── decontamination.py          Parallel sheets + orphan fragments
│   ├── detection_hypocotyle.py     Bounded collar + hypocotyl detection
│   └── root_traits_full.py         Full trait set
├── docs/traits.md              Trait reference
├── docs/limitations.md         Known limitations
├── params_example.txt          Config for the bundled example
├── example/roots/sample_S1.npy Synthetic example dataset (versioned)
├── tests/test_invariants.py    Invariant tests (pytest)
├── validation/validate_phantoms.py  Accuracy check on known-geometry phantoms
├── params.txt                  All settings
├── pyproject.toml              Package metadata (pip install -e .)
├── requirements.txt
├── makeDir.sh                  Creates data/ and results/ folders
├── data/                       Input volumes (not versioned)
└── results/                    Output (Excel, figures, checkpoint)
```

Run `run_pipeline.py` from the root folder (or `pip install -e .` first) so that `rootctrait` is importable.

## Known limitations

RootCTrait is a research tool with clearly stated limits: it depends entirely on
the upstream segmentation, oblique or poorly segmented hypocotyls can be missed,
thresholds are calibrated on soybean CT at ~0.39 x 0.39 x 0.2 mm, and there is no
ground truth for 3D roots (validation is by consistency, reproducibility and
visual inspection, not absolute accuracy). See [`docs/limitations.md`](docs/limitations.md)
for the full discussion.

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

To make the software formally citable with a permanent DOI, create a tagged
release and archive it on Zenodo:

1. Sign in to Zenodo with your GitHub account and enable archiving for this
   repository (Zenodo settings, GitHub tab).
2. On GitHub, create a release (e.g. tag `v1.0.0`). Zenodo archives it and mints
   a DOI automatically.
3. Add the DOI to `CITATION.cff` (see the commented `identifiers` block) and to
   this README, then cite it in your GWAS papers as, for example:
   "Root traits were extracted with RootCTrait v1.0 (DOI:10.5281/zenodo.XXXXXXX)."

## License

Released under the MIT License. See [`LICENSE`](LICENSE). Copyright (c) 2026
Marc-Antoine Chiasson.
