# RootCTrait trait reference

One row per sample is written to the output Excel table. Lengths are in cm,
diameters in mm, volumes in cm3, angles in degrees, ratios and counts are
dimensionless. Depth is measured from the raised collar.

## Bookkeeping (not biological traits)

| Column     | Meaning                                                        |
|------------|----------------------------------------------------------------|
| `n_brut`   | number of skeleton segments before cleaning                    |
| `n_retire` | number of segments removed as pollution                        |
| `%retire`  | fraction removed (a scan-quality indicator; use as covariate)  |

## Length

| Trait  | Definition                                              | Unit |
|--------|---------------------------------------------------------|------|
| `LRP`  | primary root (pivot) length                             | cm   |
| `TRL`  | total root length (sum of all segments)                 | cm   |
| `LTRL` | total lateral length (order >= 2)                       | cm   |
| `MLRL` | mean lateral length                                     | cm   |

## Counts and topology

| Trait            | Definition                                          | Unit  |
|------------------|-----------------------------------------------------|-------|
| `NRL`            | number of lateral roots                             | count |
| `NRL_court_<5`   | laterals shorter than 5 mm                          | count |
| `NRL_moyen_5_15` | laterals 5 to 15 mm                                 | count |
| `NRL_long_>15`   | laterals longer than 15 mm                          | count |
| `NT`             | number of tips (apices, leaf segments)              | count |
| `NBP`            | number of branch sites carrying a lateral           | count |
| `MaxO`           | maximum branching order                             | count |
| `NTR`            | number of roots emerging at the collar (basal roots)| count |

## Depth

| Trait   | Definition                                              | Unit |
|---------|---------------------------------------------------------|------|
| `PM`    | maximum rooting depth                                   | cm   |
| `D50`   | median depth of the skeleton point cloud                | cm   |
| `D95`   | 95th percentile depth                                   | cm   |
| `DMAX`  | maximum local diameter along the system                 | mm   |
| `DD_cv` | coefficient of variation of segment diameters           | ---  |

## Width and shape of the silhouette

| Trait  | Definition                                              | Unit |
|--------|---------------------------------------------------------|------|
| `WX`   | maximum width along X                                   | cm   |
| `WZ`   | maximum width along Z                                   | cm   |
| `LM`   | maximum lateral spread, max(WX, WZ)                     | cm   |
| `RLP`  | width-to-depth ratio, LM / PM                          | ---  |
| `W25`  | width at 25% of maximum depth                          | cm   |
| `W50`  | width at 50% of maximum depth                          | cm   |
| `W75`  | width at 75% of maximum depth                          | cm   |

## Angles

| Trait        | Definition                                                    | Unit |
|--------------|---------------------------------------------------------------|------|
| `ANGsys`     | system angle: length-weighted mean segment angle to vertical  | deg  |
| `ACRL`       | mean lateral angle to vertical                                 | deg  |
| `ANGO2`      | mean angle of order-2 laterals (start-to-end vector) to vertical | deg |
| `ANGO2_sd`   | standard deviation of the order-2 angles (framework regularity) | deg |
| `ANGO2_init` | mean initial angle (first 4 mm) of order-2 laterals to vertical | deg |

0 degrees = vertical, 90 degrees = horizontal.

## Volume and surface

These are computed on the raw mask and are the most sensitive to segmentation
noise and to water content in the pot.

| Trait | Definition                                              | Unit    |
|-------|---------------------------------------------------------|---------|
| `CHV` | convex hull volume of the skeleton point cloud          | cm3     |
| `VRT` | root volume (from the mask)                             | cm3     |
| `SRT` | root surface area (from the mask)                       | cm2     |
| `IC`  | compactness, VRT / CHV (<= 1)                           | ---     |
| `SRL` | specific root length, TRL / VRT                         | mm/mm3  |

## Density, spacing, diameter, tortuosity

| Trait   | Definition                                              | Unit  |
|---------|---------------------------------------------------------|-------|
| `DR`    | root density, laterals per cm of pivot                  | nb/cm |
| `IBD`   | mean inter-branch distance along the pivot              | cm    |
| `DRP`   | mean primary root diameter                              | mm    |
| `DRS`   | mean lateral diameter                                   | mm    |
| `TAPER` | pivot taper: proximal vs distal diameter, relative, per cm | 1/cm |
| `TOR`   | mean lateral tortuosity (length / straight-line distance, >= 1) | --- |

## Notes for GWAS

- Size and count traits (`TRL`, `LTRL`, `NRL`, `NT`, `NBP`, `VRT`, ...) are
  strongly correlated: they mostly describe one "system size" axis. Prefer a small
  non-redundant subset over the full list.
- Traits that aggregate over all segments (lengths, counts, `MaxO`, `DR`) are the
  most affected by residual pollution; targeted-geometry traits (`LRP`, angles) are
  more robust. Volume traits (`CHV`, `VRT`, `SRT`, `IC`, `SRL`) are the most
  sensitive.
- Include the batch and `%retire` as covariates to absorb the scan-quality
  confound. Treat `%retire`, `MaxO` and collar-related quantities as quality
  indicators rather than biological traits.
