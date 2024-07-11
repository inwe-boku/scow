# %% imports
import pandas as pd
from scow.utils import subprocess_rscript, postprocess_spatialdc, leave_one_out
from scow.tex_convert import summary_statistics_to_latex, coeffs_to_latex, stats_to_latex, convert_names
from config import repo, data_ver, num_runs, selection_criterion

compute = True

# %% generate results folder if not existent
(repo / "data" / "results").mkdir(parents=True, exist_ok=True)

# %% summary stats
df = pd.read_csv(repo / "data" / "processed" / f"dc_data_{data_ver}.csv")
sens_name = "" if data_ver == "2014" else "_modern"
summary_statistics_to_latex(df, f"sum_stats{sens_name}.tex")

# %% set model variables
binary_variables = [
    "airports_buff",
    "alps_convention",
    "broadleaved",
    "coniferous",
    "important_bird_areas",
    "protected_areas:important_bird_areas",
    "pastures",
    "preservation",
    "protected_areas",
    "restricted_military_areas",
    "water_bodies",
]
count_binary = str(len(binary_variables))

float_variables = [
    "distance_buildings_in_greenland",
    "distance_existing_turbines",
    "distance_greenland_zonings",
    "distance_other_building_land",
    "distance_power_lines",
    "distance_residential_buildings",
    "distance_roads",
    "elevation",
    "min_lcoe",
    "overnight_stays",
    "slope",
    "tree_cover_density",
    "tree_cover_density:broadleaved",
    "tree_cover_density:coniferous",
]
variables = binary_variables + float_variables


# %% functions
def estimate_spatial_dc(work_dir, data_file, num_runs, run_name, count_binaries, variables, selection_criterion="BIC",
                        compute=True):
    if compute:
        subprocess_rscript(work_dir, data_file, num_runs, run_name, count_binaries, variables)
    coeff, loglik = postprocess_spatialdc(work_dir, data_file, run_name)
    insignificant = list(coeff[coeff["pval"] > 0.1].index)
    criterion = loglik[selection_criterion]
    return coeff, criterion, insignificant, loglik


def eliminate_variable(variables_list, eliminate, integer_variables):
    surviving_integers = [ivar for ivar in variables_list if ivar in integer_variables and ivar not in eliminate]
    surviving_floats = [var for var in variables_list if var not in eliminate and var not in integer_variables]
    int_count = str(len(surviving_integers))
    return surviving_integers, surviving_floats, int_count


# %% initial estimation
coeff_init, crit_init, insignificant, loglik_init = estimate_spatial_dc(
    repo,
    data_ver,
    num_runs,
    "full",
    count_binary,
    variables,
    selection_criterion,
    compute=compute,
)

vars_int_sig, vars_float_sig, count_int_significant = eliminate_variable(variables, insignificant, binary_variables)
vars_significant = vars_int_sig + vars_float_sig

# %% initial without insignificant
coeff_sig, crit_sig, insig_sig, loglik_sig = estimate_spatial_dc(
    str(repo),
    data_ver,
    num_runs,
    "full_significant",
    count_int_significant,
    vars_significant,
    selection_criterion,
    compute=compute,
)

# %% general to specific elimination of variables
gts_iteration = 1
criterion = crit_sig
exclude = insignificant.copy()

while True:
    run_name = f"exclude_{gts_iteration}"
    vars_int, vars_float, _ = eliminate_variable(vars_significant, exclude, vars_int_sig)
    select, coeff = leave_one_out(vars_int + vars_float, vars_int, run_name, repo, data_ver, num_runs, compute=compute)

    # Find the variable to exclude that minimize the selection criterion
    best_variable_to_exclude = select.loc[select[selection_criterion].idxmin(), "variable"]
    # If excluding the best variable does not improve selection criterion, break out of the loop
    if select[selection_criterion].min() >= criterion:
        break
    # Update exclusion list and selection criterion
    exclude.append(best_variable_to_exclude)
    criterion = select[selection_criterion].min()
    # Increment iteration counter
    gts_iteration += 1

# %% read final models
full_name = "full"
full_coefs, full_stats = postprocess_spatialdc(repo, data_ver, full_name)

if gts_iteration == 1:
    parsimonious_name = "full_significant"
else:
    parsimonious_name = f"exclude_{gts_iteration - 1}_{exclude[-1].replace(':', '_')}"

with open(f"parsimonious_name{sens_name}.txt", "w") as txt_file:
    txt_file.write(parsimonious_name)

parsimonious_coeffs, parsimonious_stats = postprocess_spatialdc(repo, data_ver, parsimonious_name)

# transform to latex
latex_coeffs = coeffs_to_latex(full_coefs, parsimonious_coeffs, 5)
latex_stats = stats_to_latex(full_stats, parsimonious_stats)

# Write LaTeX code to a .tex file
with open(f"dc_estimates{sens_name}.tex", "w") as tex_file:
    tex_file.write(latex_coeffs + "\n\\midrule\n" + latex_stats)

# %% compute and write WTP to latex
ests = pd.read_csv(repo / 'data' / 'results' / f'spatialdc_coefs_{data_ver}_{parsimonious_name}.csv')
ests = ests.set_index("term")

wtp_confidence = []
for mod_nr in range(1, int(num_runs) + 1):
    wtp_confidence.append(
        ests[ests["mod"] == mod_nr]["estimate"] / ests[ests["mod"] == mod_nr].loc["min_lcoe", "estimate"])

wtp_confidence = pd.concat(wtp_confidence, axis=1)
wtp = pd.DataFrame(index=wtp_confidence.index, columns=["Estimate", "ci95_lower", "ci95_upper"])
wtp["ci95_lower"] = wtp_confidence.quantile(0.025, axis=1)
wtp["ci95_upper"] = wtp_confidence.quantile(0.975, axis=1)
wtp["Estimate"] = parsimonious_coeffs['estimate'] / parsimonious_coeffs.loc['min_lcoe', 'estimate']
wtp['95% confidence interval'] = wtp.apply(lambda row: f"[{row['ci95_lower']:.3f}, {row['ci95_upper']:.3f}]", axis=1)
wtp = wtp.drop(columns=['ci95_lower', 'ci95_upper'])


from scow.tex_convert import units
units_wtp = {
    key: r"\euro/\si{\mega\watt\hour} / " + value if value not in ["--", r"\si{\percent}"] else r"\euro/\si{\mega\watt\hour}"
    for key, value in units.items()
}
wtp.loc['coniferous:tree_cover_density', 'Unit'] = r"\euro/\si{\mega\watt\hour}"

wtp = wtp.rename(index=convert_names)
wtp['Unit'] = wtp.index.map(units_wtp)
wtp.loc['Quasi-LCOE', 'Unit'] = ""
wtp.loc['Touristic overnights', 'Unit'] = r"\euro/\si{\mega\watt\hour}" + " / 1000"
wtp.loc[wtp['Unit'].isna(), 'Unit'] = r"\euro/\si{\mega\watt\hour}"
wtp = wtp[['Unit', 'Estimate', '95% confidence interval']]
latex_wtp = wtp.to_latex(float_format="%.3f")
lines = latex_wtp.split('\n')[5:-3]  # Remove the first and last lines1
if lines[-1].endswith('\\\\'):
    lines[-1] = lines[-1][:-2]
latex_wtp_no_header_footer = '\n'.join(lines)

with open(f"valuation{sens_name}.tex", "w") as tex_file:
    tex_file.write(latex_wtp_no_header_footer)
