# %% imports
import subprocess
from config import repo

# %% data preprocessing
subprocess.run(f"python {repo}/01_spatial_data.py")
# import src.data_download

# %% estimation of willingness to pay
subprocess.run(f"python {repo}/02_estimate_wtp.py")

# %% postprocessing of estimation results
subprocess.run(f"python {repo}/03_postprocess_wtp.py")

# %% optimal choice of wind turbine sites
subprocess.run(f"python {repo}/04_optimal_siting.py")

# %% postprocessing of wind turbine sites
subprocess.run(f"python {repo}/05_postprocess_sites.py")
