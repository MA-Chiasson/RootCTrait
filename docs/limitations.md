# Known limitations

RootCTrait is an honest research tool. The points below should be reported
alongside any results.

## 1. It depends entirely on the upstream segmentation

RootCTrait works on an already segmented binary mask. It does not segment the CT
volume itself. Segmentation errors propagate directly into the traits:

- Surface "sheets" / dense pollution near the collar inflate any trait that
  aggregates over segments (total length, counts, branching order, density).
- Volume traits computed on the raw mask (`CHV`, `VRT`, `SRT`, `IC`, `SRL`) are the
  most sensitive, and are affected by water content in the pot at scan time.

The decontamination step removes a large part of this, but it cannot recover
information that the segmentation lost or invented. Use `%removed` as a
scan-quality covariate, and inspect flagged samples with the 3D figures.

## 2. Oblique or poorly segmented hypocotyls can be missed

The hypocotyl is detected as a vertical column above the collar. A hypocotyl that
leans strongly, or that is broken up by the segmentation, is geometrically hard to
distinguish from a root and may be partially kept. Small missed stubs have a
limited effect on the traits; larger oblique hypocotyls are worth a visual check.

## 3. Calibration and scope

Thresholds (decontamination, collar climb, hypocotyl angle) were tuned and
validated on soybean CT scans at a voxel size of about 0.39 x 0.39 x 0.2 mm. For a
different species, resolution, or imaging modality, revalidate visually and adjust
the parameters in `params.txt` if needed. The multi-format loader is implemented
for `.mat`, `.tif`, `.npy` and `.nii`, but has been tested most thoroughly on
`.mat`.

## 4. No ground truth for 3D roots

There is no reference measurement for the true architecture of a 3D root system.
As for every tool in this field, RootCTrait can be assessed for internal
consistency, reproducibility across replicates, and visual validation, but not for
absolute accuracy. Trait values should be read as reproducible descriptors, not as
exact ground-truth measurements.
