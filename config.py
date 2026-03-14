# %% imports
from pathlib import Path
import shutil

# %% config
repo = Path(__file__).resolve().parent
atlas_path = repo / "data" / "atlas"
data_ver = "2014"  # "2014" or "modern"

num_runs = "500"
selection_criterion = "BIC"
replication_seed = 42

spacing = 3  # retained default discrete spacing rule for the accepted Python-native stage-4 path
nturbines = "auto"

rdir = Path(shutil.which("Rscript") or "Rscript")
