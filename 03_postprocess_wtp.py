from scow.logging_utils import configure_script_logging
from scow.postprocess_wtp_helpers import (
    build_interaction_layers,
    compute_social_cost_dataset,
)
from scow.paths import (
    parsimonious_name_path,
    soco_path,
)
from scow.stages import load_atlas_for_analysis
from scow.utils import postprocess_spatialdc
from config import repo, atlas_path, data_ver

logger = configure_script_logging("stage3")

with open(parsimonious_name_path(data_ver), "r") as f:
    parsimonious_name = f.read()

# %% load estimated model
parsimonious_coeffs, _ = postprocess_spatialdc(repo, data_ver, parsimonious_name)
parsimonious_coeffs["valuation"] = parsimonious_coeffs["estimate"] / parsimonious_coeffs.loc["min_lcoe", "estimate"]

# %% load site data
atlas = load_atlas_for_analysis({"atlas_path": atlas_path, "scenario": data_ver})

# generate interaction data
landscape_data = build_interaction_layers(atlas.landscape.data, list(parsimonious_coeffs.index))

if "min_lcoe" not in landscape_data.data_vars:
    landscape_data = landscape_data.merge(atlas.wind.data["min_lcoe"].compute())
spatial_data = landscape_data[list(parsimonious_coeffs.index)]

# %% compute social cost
soco = compute_social_cost_dataset(spatial_data, parsimonious_coeffs["valuation"])
soco.to_netcdf(soco_path(data_ver))
logger.info("Wrote social-cost dataset to %s", soco_path(data_ver))
