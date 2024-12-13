# %% imports
from pathlib import Path

# %% config
repo = Path("/Users/nwesec/repos/scow")
data_ver = "2014"  # "2014" or "modern"

num_runs = "2500"
selection_criterion = "BIC"

spacing = 3
nturbines = "auto"

rdir = Path("/usr/local/bin/Rscript")
gamsdir = Path("/Library/Frameworks/GAMS.framework/Versions/46/Resources/gams")

gams_conf = {
    'gams_model': repo / "opt" / "location_selection.gms",
    'gdx_input': repo / "opt" / "input_data.gdx",
    'gdx_output': repo / "opt",
    'gams_exe': gamsdir.parent
}


