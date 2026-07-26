# Data

The individual dengue case data used in this study are **not** included in this
repository. They were obtained from collaborating hospitals under ethical
approval and data-sharing restrictions and contain identifiable patient
information; they cannot be redistributed. See the *Data availability* statement
in the paper.

To reproduce the analysis, place the following inputs in this `data/` folder:

| File | Description | Source |
|------|-------------|--------|
| `Individual dengue data.xlsx` | Case-level records (reporting hospital, subdistrict of residence, travel distance) plus the `total case by hospital` and `total case by subdistrict` summary sheets | Collaborating hospitals (restricted) |
| `Dengue Spatial Data 5d.xlsx` | Subdistrict population density (used to form the population offset) | BPS "Makassar in Figures" |
| `shapefile/Makassar_no_sangkarrang.shp` (+ `.shx`, `.dbf`, `.prj`) | Administrative boundaries of the 14 mainland subdistricts | Public administrative boundary data |

The scripts read from `data/` via the `DATA_DIR` variable set in `R/00_setup.R`.
Only the two aggregated summary tables (`total case by subdistrict`,
`total case by hospital`) and the offset are strictly required to run the
area-level models; the case-level sheet is required for the full point-process
model and the origin-recovery cross-validation.
