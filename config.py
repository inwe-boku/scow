# %% imports
from pathlib import Path

# %% config
repo = Path("c:/git_repos/repscow")
data_ver = "2014"  # "2014" or "modern"

num_runs = "2500"
selection_criterion = "BIC"

spacing = 3
nturbines = "auto"

rdir = Path("C:/myprogs/R/R-4.2.2/bin/Rscript")
gamsdir = Path("c:/myprogs/GAMS/45")

gams_conf = {
    'gams_model': repo / "opt" / "location_selection.gms",
    'gdx_input': repo / "opt" / "input_data.gdx",
    'gdx_output': repo / "opt",
    'gams_exe': gamsdir
}


