# Inferring the Local Social Cost of Wind Turbines

This project estimates the valuation of spatial attributes as implied by authorities' renewable-energy zoning decisions. Data processing uses the [cleo package](https://github.com/sebwehrle/cleo).

## Study Coverage

The repository now tracks three explicit case-study targets:

- `lower_austria_wind_2014`
- `lower_austria_pv_2022`
- `lower_austria_wind_2024`

Current implementation status:

- only `lower_austria_wind_2014` is implemented end-to-end in the numbered pipeline
- `lower_austria_wind_2024` now has an explicit live stage-1 payload and explicit stage-2 to stage-5 artifact namespaces, but it is still not marked executable because the 2024 covariate recipe and missing-input inventory are not yet closed out
- `lower_austria_pv_2022` has core public documents ingested, but no dedicated stage-1 data bundle or code path yet

Important clarification:

- the older runtime label `modern` is a legacy internal scenario label in parts of the current codebase
- it is not a canonical study id
- the legacy runtime modes `2014` and `modern` both remain supported for direct workflow comparisons
- explicit study ids keep a separate artifact namespace and do not accept conflicting runtime-scenario overrides
- the tracked high-level status file is [PROJECT_STATUS.md](/Users/nwesec/code/scow/PROJECT_STATUS.md)

## Data
The analysis includes proprietary data which must be obtained from:
- **Important Bird Areas**: BirdLife International via [BirdLife Data Zone](https://datazone.birdlife.org/site/requestgis)
- **Overnight Stays in Austrian Municipalities**: Statistics Austria (request via [info@statistik.gv.at](mailto:info@statistik.gv.at))

All other data are downloaded from sources referenced in:
```python
from scow.site_data import generate_data_dict
```

For the World Database on Protected Areas (WDPA), if the provided link is broken, find the latest link at [Protected Planet](https://www.protectedplanet.net/country/AUT).

## Installation
The code is developed in a Python 3.12 environment and relies heavily on the local `cleo` package, which currently requires Python 3.11 or newer. Other Python dependencies are listed in `environment.yml`. To create the canonical Python environment:
```bash
conda env create -n scow -f environment.yml
```
For `cleo` installation and usage, refer to the [cleo GitHub page](https://github.com/sebwehrle/cleo).

Apart from python, installations of GDAL, [R](https://www.r-project.org/) and a working [Gurobi](https://www.gurobi.com/) setup for the Python environment are also required.

Stage 2 now expects its R package set to be restored ahead of execution. The tracked lockfile is `renv.lock`. The canonical bootstrap path is:

```bash
Rscript --vanilla -e 'install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(lockfile="renv.lock")'
```

`discrete_choice.R` now bootstraps the project library from that lockfile at runtime, so later stage-2 runs and the focused stage-2 test use the tracked package set rather than whatever happens to be installed globally.

The repository also tracks the standard `renv` activation files:

- `.Rprofile`
- `renv/activate.R`

Those files support interactive and non-`--vanilla` R sessions. The canonical stage-2 pipeline still calls `Rscript --vanilla`, so the runtime guarantee for numbered-script execution remains the explicit bootstrap inside `discrete_choice.R`.

The current locked stage-2 package surface includes at least:
- `logitr`
- `tidyverse`
- `broom`
- `jsonlite`

To run the focused validation suite from the Python environment:

```bash
pytest
```

The suite includes guarded checks for the stage-2 R runtime and the optional Gurobi solver path. Those checks will skip rather than fail when the corresponding runtime is not provisioned locally.

## Execution
### Configuration
Configure the following settings in `config.py`:
```python
from pathlib import Path
import shutil

from scow.study import ScowStudy

repo = Path(__file__).resolve().parent
study = ScowStudy.from_study_id(
    "lower_austria_wind_2014",
    num_runs=500,              # repeated discrete-choice iterations
    selection_criterion="BIC", # "AIC" or "BIC"
    seed=42,                   # stage-2 replication seed
    spacing=3,                 # stage-4 spacing rule in pixels
    num_turbines="auto",       # stage-4 turbine-count rule
)

data_ver = study.scenario
artifact_key = study.artifact_key
num_runs = str(study.num_runs)
selection_criterion = study.selection_criterion
replication_seed = study.seed
spacing = study.spacing
nturbines = study.num_turbines

rdir = Path(shutil.which("Rscript") or "Rscript")  # resolved from PATH at runtime
```

The selected `study` is now the authoritative home for the current project defaults and assumptions. The legacy module-level aliases such as `data_ver`, `num_runs`, `spacing`, and `nturbines` are still exported for compatibility with the current scripts.

Contract note:

- `study.scenario` controls runtime semantics such as turbine family and atlas loading
- `study.artifact_key` controls the artifact namespace written by stages 2 to 5
- in legacy scenario mode, `artifact_key == scenario`
- in explicit study mode, `artifact_key` is the study id while `scenario` remains the runtime label

For explicit study-level project planning rather than the legacy runtime scenario switch, use:

- `scow.study_registry.get_study(...)`
- `scow.pipeline.build_pipeline(...)`

Those helpers expose the three tracked study ids listed above and clarify which study is executable today.

Example:

```python
from scow.pipeline import build_pipeline

pipeline = build_pipeline("lower_austria_wind_2024")
print(pipeline.executable)      # False
print(pipeline.stages[0].status)  # "implemented"
```

To run the live 2024 stage-1 payload only:

```python
from scow.study import ScowStudy

study = ScowStudy.from_study_id("lower_austria_wind_2024")
```

`data/results` and `figures` are created automatically by the path helpers when downstream stages write outputs.

### Data Processing
Execute the scripts sequentially:

1. **01_spatial_data.py**: Retrieves and preprocesses wind resources data from the Global Wind Atlas and other spatial characteristics. The `lower_austria_wind_2014` runtime still writes `dc_data_2014.csv` and `dc_data_modern.csv` in `data/processed`. The explicit `lower_austria_wind_2024` stage-1 payload writes `dc_data_lower_austria_wind_2024.csv`.

2. **02_estimate_wtp.py**: A Python-first stage-2 driver around R code in `discrete_choice.R`, producing files like `spatialdc_coefs_{datafile}_{run_name}.csv` in `data/results`.

   Stage 2 now also writes `spatialdc_meta_{datafile}_{run_name}.json` sidecars so repeated model-search runs only reuse cached artifacts when the artifact key, seed, run count, integer partition, and exact term list still match.

   The primary configuration surface is now the Python `EstimationSpec` / preset API in `scow.estimation`. The script still runs directly with no arguments, and also supports a thin CLI for common overrides, for example:

   ```bash
   python 02_estimate_wtp.py --preset baseline --single --num-runs 50
   ```

   For explicit studies, conflicting `--scenario` overrides are rejected. For legacy scenario mode, `2014` and `modern` remain valid parallel runs.

   For noncanonical experimental runs from Python, use an explicit `output_label` so experimental artifacts do not overwrite the canonical baseline namespace.

   The legacy `pval` / significance-star outputs are retained for compatibility with the current tables and selection loop, but they summarize coefficient stability across repeated negative samples rather than classical single-fit logit inference.

3. **03_postprocess_wtp.py**: Computes the local social-cost artifact and writes `soco.nc` in `data/results`.

4. **04_optimal_siting.py**: Solves the location optimization problem through the Python-native `gurobipy` stage-4 path, reusing validated chunk caches under `opt/` when possible and writing final results to `opt_locations_{artifact_key}_{objective}.csv` and `opt_cost_{artifact_key}_{objective}.csv` in `data/results`.

   Stage 4 now also writes `opt_meta_{artifact_key}_{objective}.json` sidecars and uses shared validators for the final `opt_locations_*` / `opt_cost_*` / `opt_meta_*` bundle contract.

5. **05_postprocess_sites.py**: Unified reporting stage for both WTP/social-cost outputs and site-allocation outputs. It reads `soco_{artifact_key}.nc` plus the matching stage-4 outputs, then writes the WTP maps/distributions, turbine-allocation maps, social-cost curve, and LaTeX cost comparison table under the same artifact key. The 2014 baseline still keeps the canonical unsuffixed root outputs.

## Status Tracking

The repository now keeps the high-level study and implementation status in:

- [PROJECT_STATUS.md](/Users/nwesec/code/scow/PROJECT_STATUS.md)

The explicit study registry and pipeline-builder surface live in:

- [scow/study_registry.py](/Users/nwesec/code/scow/scow/study_registry.py)
- [scow/pipeline.py](/Users/nwesec/code/scow/scow/pipeline.py)

### Validation
The repository now includes a small focused validation suite covering:

- stage-2 reproducibility
- stage-2 artifact handoff on a tiny synthetic design matrix
- stage-3 social-cost dataset contract checks
- stage-4 final-output bundle validation and loading
- stage-4 chunk candidate/conflict semantics
- stage-5 generation-target accounting
- stage-4 to stage-5 artifact handoff on synthetic outputs
- IGW scraper cache-contract validation
- on-demand output-path creation

### Logging
The numbered pipeline scripts now share a lightweight stderr logger with timestamped stage labels. This keeps the default workflow script-driven and console-oriented while making long stage runs and artifact writes easier to trace. The repository does not enable mandatory per-run log files by default.
