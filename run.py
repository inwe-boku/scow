# %% imports
from pathlib import Path
import subprocess
import sys
from config import repo
from scow.logging_utils import configure_script_logging

scripts = [
    "01_spatial_data.py",
    "02_estimate_wtp.py",
    "03_postprocess_wtp.py",
    "04_optimal_siting.py",
    "05_postprocess_sites.py",
]

logger = configure_script_logging("run")

for script in scripts:
    logger.info("Launching %s", script)
    subprocess.run([sys.executable, str(Path(repo) / script)], check=True)
    logger.info("Completed %s", script)
