# %% imports
import subprocess
import sys
from config import repo

scripts = [
    "01_spatial_data.py",
    "02_estimate_wtp.py",
    "03_postprocess_wtp.py",
    "04_optimal_siting.py",
    "05_postprocess_sites.py",
]

for script in scripts:
    subprocess.run([sys.executable, str(repo / script)], check=True)
