# RootCTrait trait reference

One row per sample is written to the output Excel table: 38 traits plus three
bookkeeping columns. Lengths are in cm, diameters in mm, volumes in cm3, angles in
degrees; ratios and counts are dimensionless. Unless stated otherwise, traits are
computed on the cleaned root system (hypocotyl, sheet artifacts and detached
fragments excluded), and depths are measured from the raised collar.

Lateral counts and lengths are computed on **segments**: a segment is a maximal
chain of skeleton nodes of the same branching order, so a branched lateral
contributes one segment per order.

## Bookkeeping (not biological traits)

| Column      | Meaning                                                        |
|-------------|----------------------------------------------------------------|
| `n_raw`     | number of skeleton segments before cleaning                    |
| `n_removed` | number of segments removed by decontamination                  |
| `%removed`  | fraction removed (a scan-quality indicator; use as covariate)  |
| `pivot_return` | upward return of the pivot after its deepest point, mm (quality control; check samples above 3 mm) |

## Length

| Trait  | Definition                                                                  | Unit |
|--------|-----------------------------------------------------------------------------|------|
| `LRP`  | primary root (pivot) length, along the whole pivot path from the raised collar to its tip | cm |
| `TRL`  | total root length (sum of all segment lengths)                              | cm   |
| `LTRL` | total lateral length (segments of order >= 2)                               | cm   |
| `MLRL` | mean length of lateral segments (order >= 2)                                | cm   |

## Counts and topology

| Trait             | Definition                                                         | Unit  |
|-------------------|--------------------------------------------------------------------|-------|
| `NRL`             | number of lateral segments (order >= 2)                            | count |
| `NRL_short_<5`    | lateral segments shorter than 5 mm                                 | count |
| `NRL_medium_5_15` | lateral segments 5 to 15 mm                                        | count |
| `NRL_long_>15`    | lateral segments 15 mm or longer                                   | count |
| `NT`              | number of terminal segments (apices)                               | count |
| `NBP`             | number of segments bearing at least one lateral                    | count |
| `MaxO`            | maximum branching order                                            | count |
| `NTR`             | proxy count of roots emerging near the collar: segments of order <= 2 with an end within 4 mm of the raised collar (minimum 1) | count |
| `DR`              | lateral density: `NRL` per cm of `LRP`                             | nb/cm |
| `IBD`             | mean distance along the pivot between successive lateral insertion points | cm |

## Depth

Skeleton points more than 2 mm above the raised collar are ignored.

| Trait | Definition                                                   | Unit |
|-------|--------------------------------------------------------------|------|
| `PM`  | maximum rooting depth below the raised collar                | cm   |
| `D50` | median depth of the skeleton point cloud                     | cm   |
| `D95` | 95th percentile depth of the skeleton point cloud            | cm   |

## Width and shape of the silhouette

| Trait | Definition                                                                  | Unit |
|-------|-----------------------------------------------------------------------------|------|
| `WX`  | maximum extent of the skeleton along X                                      | cm   |
| `WZ`  | maximum extent of the skeleton along Z                                      | cm   |
| `LM`  | maximum lateral spread, max(`WX`, `WZ`)                                     | cm   |
| `W25` | width at 25% of `PM`: larger of the X and Z extents within a ± 2 mm depth band | cm |
| `W50` | width at 50% of `PM`, as `W25`                                              | cm   |
| `W75` | width at 75% of `PM`, as `W25`                                              | cm   |
| `RLP` | width to depth ratio, `LM` / `PM`                                           | ---  |

## Angles

0 degrees = vertical, 90 degrees = horizontal.

| Trait      | Definition                                                                  | Unit |
|------------|-----------------------------------------------------------------------------|------|
| `ANGO2`    | mean angle to the vertical of order 2 segments (end to end vector)          | deg  |
| `ANGO2_sd` | standard deviation of the order 2 angles to the vertical                    | deg  |
| `ANGI`     | mean insertion angle of laterals (order >= 2) relative to their parent: angle between the chord of the first 4 mm of the lateral and the local axis of the parent at the attachment point, folded into 0 to 90 | deg |

Three earlier angle traits (`ANGsys`, `ACRL`, `ANGO2_init`) were dropped as
redundant with `ANGO2` and are no longer computed.

## Volume and surface

`VRT` and `SRT` are computed on the part of the binarized mask attached to the
cleaned roots: each mask voxel is assigned to its nearest skeleton voxel, and kept
when that voxel belongs to the cleaned root system. Mask components that are not
connected to the cleaned roots are excluded first. The hypocotyl, removed layers,
orphan and isolated fragments are therefore excluded, consistently with the
skeleton traits.
They remain the most sensitive traits to segmentation and to water content in the
pot.

| Trait | Definition                                                         | Unit   |
|-------|--------------------------------------------------------------------|--------|
| `CHV` | convex hull volume of the cleaned skeleton point cloud             | cm3    |
| `VRT` | root volume: voxel count of the cleaned root mask times the voxel volume | cm3 |
| `SRT` | root surface area: area of the marching cubes mesh of the cleaned root mask | cm2 |
| `IC`  | compactness, `VRT` / `CHV` (set to missing when above 1)           | ---    |
| `SRL` | length per root volume, `TRL` / `VRT`                              | cm/cm3 |

## Diameter and form

| Trait   | Definition                                                                  | Unit |
|---------|-----------------------------------------------------------------------------|------|
| `DRP`   | mean primary root diameter, mean of 2 edt along the corrected pivot         | mm   |
| `DRS`   | mean lateral diameter, mean over lateral segments of 2 median(edt)          | mm   |
| `DMAX`  | maximum local diameter, 2 max(edt) over the cleaned skeleton                | mm   |
| `DD_cv` | coefficient of variation of segment diameters (all orders)                  | ---  |
| `TAPER` | pivot taper: (proximal 25% diameter − distal 25% diameter) / proximal diameter, per cm of pivot | 1/cm |
| `TOR`   | mean lateral tortuosity (segment length / straight line distance, >= 1)     | ---  |

## Notes for GWAS

- Size and count traits (`TRL`, `LTRL`, `NRL`, `NT`, `NBP`, `VRT`, ...) are
  strongly correlated: they mostly describe one "system size" axis. Prefer a small
  non-redundant subset over the full list.
- Traits that aggregate over all segments (lengths, counts, `MaxO`, `DR`) are the
  most affected by residual pollution; targeted-geometry traits (`LRP`, angles) are
  more robust. Mask based traits (`VRT`, `SRT`, `IC`, `SRL`, and to a lesser extent
  `CHV`) are the most sensitive.
- Include the batch and `%removed` as covariates to absorb the scan-quality
  confound. Treat `%removed`, `MaxO` and collar-related quantities as quality
  indicators rather than biological traits.
