# %% imports
from pathlib import Path
import shutil

from scow.study import ScowStudy

# %% config
repo = Path(__file__).resolve().parent
atlas_path = repo / "data" / "atlas"
study = ScowStudy.from_study_id("lower_austria_wind_2014")
data_ver = study.scenario  # legacy runtime scenario label used by the current executable workflow

num_runs = str(study.num_runs)
selection_criterion = study.selection_criterion
replication_seed = study.seed

spacing = study.spacing  # retained default discrete spacing rule for the accepted Python-native stage-4 path
nturbines = study.num_turbines

rdir = Path(shutil.which("Rscript") or "Rscript")
