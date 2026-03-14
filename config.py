# %% imports
from pathlib import Path
import shutil

from scow.study import ScowStudy

# %% config
repo = Path(__file__).resolve().parent
atlas_path = repo / "data" / "atlas"
study = ScowStudy(scenario="2014")
data_ver = study.scenario  # "2014" or "modern"

num_runs = str(study.num_runs)
selection_criterion = study.selection_criterion
replication_seed = study.seed

spacing = study.spacing  # retained default discrete spacing rule for the accepted Python-native stage-4 path
nturbines = study.num_turbines

rdir = Path(shutil.which("Rscript") or "Rscript")
