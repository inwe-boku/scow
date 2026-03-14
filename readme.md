# Inferring the Local Social Cost of Wind Turbines

This project estimates the valuation of spatial attributes as implied by the Lower Austrian wind power zoning. Data processing uses the [cleo package](https://github.com/sebwehrle/cleo).

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

repo = Path(__file__).resolve().parent
data_ver = "2014"  # "2014" for the time of zoning decision or "modern" for current wind turbines

num_runs = "500"  # number of iterations for discrete choice estimation
selection_criterion = "BIC"  # options: "AIC" or "BIC"
replication_seed = 42  # seed for repeated discrete-choice estimation

spacing = 3  # minimum pixels between wind turbines for optimal location
nturbines = "auto"  # automatic selection of deployed wind turbines

rdir = Path(shutil.which("Rscript") or "Rscript")  # resolved from PATH at runtime
```

`data/results` and `figures` are created automatically by the path helpers when downstream stages write outputs.

### Data Processing
Execute the scripts sequentially:

1. **01_spatial_data.py**: Retrieves and preprocesses wind resources data from the Global Wind Atlas and other spatial characteristics, saving `dc_data_2014.csv` and `dc_data_modern.csv` in `data/processed`.

2. **02_estimate_wtp.py**: A Python-first stage-2 driver around R code in `discrete_choice.R`, producing files like `spatialdc_coefs_{datafile}_{run_name}.csv` in `data/results`.

   The primary configuration surface is now the Python `EstimationSpec` / preset API in `scow.estimation`. The script still runs directly with no arguments, and also supports a thin CLI for common overrides, for example:

   ```bash
   python 02_estimate_wtp.py --preset baseline --single --num-runs 50
   ```

   For noncanonical experimental runs from Python, use an explicit `output_label` so experimental artifacts do not overwrite the canonical baseline namespace.

3. **03_postprocess_wtp.py**: Computes the local social-cost artifact and writes `soco.nc` in `data/results`.

4. **04_optimal_siting.py**: Solves the location optimization problem through the Python-native `gurobipy` stage-4 path, reusing validated chunk caches under `opt/` when possible and writing final results to `opt_locations_{data_ver}_{objective}.csv` and `opt_cost_{data_ver}_{objective}.csv` in `data/results`.

   Stage 4 now also writes `opt_meta_{data_ver}_{objective}.json` sidecars and uses shared validators for the final `opt_locations_*` / `opt_cost_*` / `opt_meta_*` bundle contract.

5. **05_postprocess_sites.py**: Unified reporting stage for both WTP/social-cost outputs and site-allocation outputs. It reads `soco.nc` and stage-4 cost outputs, then writes the WTP maps/distributions, turbine-allocation maps, social-cost curve, and LaTeX cost comparison table.

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
