# vala
Estimating the value of landscape characteristics implied by wind turbine zoning in Lower Austria

Create a virtual environment ``conda env create -f environment.yml`` and activate the new environment.

In addition, the following packages need to be installed manually:
* gamstransfer via ``pip install gamsapi[transfer]==xx.y.z`` where ``xx.y.z`` represents your installed GAMS version 
  number (e.g., 45.5.0)

After package installation, ``config.py`` must be updated to point to the correct repo-directory
```
ROOTDIR = Path('c:/git_repos/scow')
GAMSDIR = Path('c:/myprogs/GAMS/45')
RDIR = Path('c:/myprogs/R/R-4.2.2')
```

# TODO: add 'data/power_curves/powercurves_research.csv' to repo

# Required non-public datasets:
* Annual tourstic overnights per municipality: 'tourism/Tabelle 30 GEH 2018.xlsx' from 
* Important bird areas from Bird Life: 'c:\git_repos\scow\data\iba\EUROPE_KBA.shp'
