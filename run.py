# %% imports
import subprocess
from config import ROOTDIR, RDIR

# %% data download
subprocess.run(f"python {ROOTDIR}/src/data_download.py")
# import src.data_download

# %% preprocessing of global wind atlas and other data
subprocess.run(f"python {ROOTDIR}/src/preprocessing.py")
# import src.preprocessing
# requires 32 GB RAM

# %% calculate capacity factors
subprocess.run(f"python {ROOTDIR}/src/process_gwa.py")
# implements bias correction
# implements investment-cost scaling for global cost only

# %% preprocess geodata for estimation
subprocess.run(f"python {ROOTDIR}/src/process_geodata.py")

# %% do estimation in R
subprocess.run(f"""{str(RDIR / "bin/Rscript.exe")} {str(ROOTDIR / "src/logitr_parallel.R")} {str(ROOTDIR)}""")
# subprocess.run(f"{RDIR}/bin/Rscript {ROOTDIR}/src/mc_choice_estimation.R {ROOTDIR}")

# %% postprocessing of R estimates


# %% location optimizations
subprocess.run(f"python {ROOTDIR}/src/optimal_locations.py")
# import src.optimal_locations


# %% plotting
subprocess.run(f"python {ROOTDIR}/src/plotting.py")
# import src.plotting

# %% compile LaTeX
subprocess.run("pdflatex xxx.tex")
