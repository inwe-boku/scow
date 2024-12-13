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
The code is developed using `python 3.9` and relies heavily on the `cleo` package. Other dependencies are listed in `environment.yml`. To create a Python environment named `myenv` with these dependencies:
```bash
conda env create -n myenv -f environment.yml
```
For `cleo` installation and usage, refer to the [cleo GitHub page](https://github.com/sebwehrle/cleo).

Apart from python, installations of GDAL, [R](https://www.r-project.org/)  and [GAMS](https://www.gams.com/) are also required.

## Execution
### Configuration
Configure the following settings in `config.py`:
```python
repo = Path("path/to/repository")
data_ver = "2014"  # "2014" for the time of zoning decision or "modern" for current wind turbines

num_runs = "2500"  # number of iterations for discrete choice estimation
selection_criterion = "BIC"  # options: "AIC" or "BIC"

spacing = 3  # minimum pixels between wind turbines for optimal location
nturbines = "auto"  # automatic selection of deployed wind turbines

rdir = Path("C:/myprogs/R/R-4.2.2/bin/Rscript")  #  path to Rscript.exe (Windows)
gamsdir = Path("c:/myprogs/GAMS/45")  # path to folder containing gams.exe (Windows)
```

### Data Processing
Execute the scripts sequentially:

1. **01_spatial_data.py**: Retrieves and preprocesses wind resources data from the Global Wind Atlas and other spatial characteristics, saving `dc_data_2014.csv` and `dc_data_modern.csv` in `data/processed`.

2. **02_estimate_wtp.py**: A Python wrapper around R code in `discrete_choice.R`, producing files like `spatialdc_coefs_{datafile}_{run_name}.csv` in `data/results`.

3. **03_postprocess_wtp.py**: Computes social costs, plots maps of private and social costs, and the distribution of various costs.

4. **04_optimal_siting.py**: Solves the location optimization problem, saving results to `opt_locations_{data_ver}_{objective}.csv` and `opt_cost_{data_ver}_{objective}.csv` in `data/results`.

5. **05_postprocess_sites.py**: Processes optimal wind turbine sites, computes local and total social costs, and plots maps of optimal locations and the social cost curve.

### Logging
Logs are kept in `data/logfile.log`.