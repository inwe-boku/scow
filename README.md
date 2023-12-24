# scow
Estimating the social cost of wind power from zoning descision in Lower Austria

Create a virtual environment ``conda env create -f environment.yml`` and activate the new environment.

In addition, the following packages need to be installed manually:
* gamstransfer via ``pip install gamsapi[transfer]==xx.y.z`` where ``xx.y.z`` represents your installed GAMS version 
  number (e.g., 45.5.0)

After package installation, ``config.py`` must be updated to point to the correct repo-directory, e.g.
```
ROOTDIR = Path('c:/git_repos/scow')
GAMSDIR = Path('c:/myprogs/GAMS/45')
RDIR = Path('c:/myprogs/R/R-4.2.2')
```

# Required non-public datasets:
* Annual tourstic overnights per municipality from Statistik Austria
* Important bird areas from Bird Life Austria / Bird Life International
