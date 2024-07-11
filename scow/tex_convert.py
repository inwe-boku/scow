# %%  imports
import pandas as pd

# %%
convert_names = {
    'zoning': 'Wind power zoning',
    'airports_buff': 'Airports',
    'alps_convention': 'Alpine Convention',
    'protected_areas': 'Protected areas',
    'important_bird_areas': 'Important bird areas',
    'important_bird_areas:protected_areas': 'Protected areas $\\times$ Important bird areas',
    'preservation': 'Landscapes worth preserving',
    'pastures': 'Pastures',
    'restricted_military_areas': 'Restricted military areas',
    'water_bodies': 'Water bodies',
    'broadleaved': 'Dom. Leaf Type: Broadleaved',
    'coniferous': 'Dom. Leaf Type: Coniferous',
    'tree_cover_density': 'Tree cover density',
    'distance_buildings_in_greenland': 'd(Greenland build. w. pres.)',
    'distance_greenland_zonings': 'd(Greenland zoning)',
    'distance_power_lines': 'd(Grid infrastructure)',
    'distance_roads': 'd(High-level roads)',
    'distance_other_building_land': 'd(Other building land)',
    'distance_existing_turbines': 'd(Pre-existing wind turbines)',
    'distance_residential_buildings': 'd(Residential buildings)',
    'elevation': 'Elevation',
    'slope': 'Terrain slope',
    'overnight_stays': 'Touristic overnights',
    'min_lcoe': 'Quasi-LCOE',
    'broadleaved:tree_cover_density': 'Tree cover density $\\times$ Broadleaved',
    'coniferous:tree_cover_density': 'Tree cover density $\\times$ Coniferous',
    }

units = {
        'Wind power zoning': '--',
        'Airports': '--',
        'Alpine Convention': '--',
        'Dom. Leaf Type: Broadleaved': '--',
        'Dom. Leaf Type: Coniferous': '--',
        'd(Greenland build. w. pres.)': r'\si{\kilo\metre}',
        'd(Pre-existing wind turbines)': r'\si{\kilo\metre}',
        'd(Greenland zoning)': r'\si{\kilo\metre}',
        'd(Other building land)': r'\si{\kilo\metre}',
        'd(Grid infrastructure)': r'\si{\kilo\metre}',
        'd(Residential buildings)': r'\si{\kilo\metre}',
        'd(High-level roads)': r'\si{\kilo\metre}',
        'Elevation': r'\si{\metre}',
        'Important bird areas': '--',
        'Quasi-LCOE': r'\euro/\si{\mega\watt\hour}',
        'Touristic overnights': "'000",
        'Pastures': '--',
        'Landscapes worth preserving': '--',
        'Protected areas': '--',
        'Restricted military areas': '--',
        'Terrain slope': r'\si{\degree}',
        'Tree cover density': r'\si{\percent}',
        'Water bodies': '--',
}

col_names = {
        'index': 'Variable',
        'unit': 'Unit',
        'mean': 'Mean',
        'std': 'Std. Dev.',
        'min': 'Min',
        'max': 'Max',
}

# %% functions


def coeffs_to_latex(full_coefs, parsimonious_coefs, significant_digits=2):
    # merge both models
    model_df = pd.merge(full_coefs, parsimonious_coefs, right_index=True, left_index=True,
                        suffixes=('_full', '_parsimonious'), how='outer')

    model_df['full'] = model_df.apply(lambda row: f"$\\underset{{({row['std_full']:.{significant_digits}f})}}{{{row['estimate_full']:.{significant_digits}f}}}^{{{row['significance_full']}}}$" if not pd.isnull(row.name) else '', axis=1)
    model_df['parsimonious'] = model_df.apply(lambda row: f"$\\underset{{({row['std_parsimonious']:.{significant_digits}f})}}{{{row['estimate_parsimonious']:.{significant_digits}f}}}^{{{row['significance_parsimonious']}}}$" if not pd.isnull(row.name) else '', axis=1)

    model_df.drop(['estimate_full', 'std_full', 'tstat_full', 'pval_full', 'significance_full', 'valuation_full',
                   'estimate_parsimonious', 'std_parsimonious', 'tstat_parsimonious', 'pval_parsimonious',
                   'significance_parsimonious', 'valuation_parsimonious'], axis=1, inplace=True)

    model_df.replace(r"$\underset{(nan)}{nan}^{nan}$", '--', inplace=True)
    model_df.rename(index=convert_names, inplace=True)

    latex_table = model_df.to_latex(escape=False)
    lines = latex_table.split('\n')[5:-3]  # Remove the first and last lines
    if lines[-1].endswith('\\\\'):
        lines[-1] = lines[-1][:-2]
    latex_table_without_header_footer = '\n'.join(lines)

    return latex_table_without_header_footer


def stats_to_latex(full_stats, parsimonious_stats):
    """

    :param full_stats:
    :param parsimonious_stats:
    :return:
    """
    def custom_float_format(x):
        if abs(x) >= 10:
            return "%.1f" % x
        else:
            return "%.3f" % x

    stats = pd.concat([full_stats, parsimonious_stats], axis=1, keys=["full", "parsimonious"])
    stats = stats.drop(["mod", "null.logLik", "AIC", "r.squared", "nobs"])
    stats = stats.rename(index={"logLik": "Log-Likelihood", "adj.r.squared": "Adj. R$^2$"})
    latex_stats = stats.to_latex(float_format=custom_float_format)
    stat_lines = latex_stats.split('\n')[4:-3]  # Remove the first and last lines
    if stat_lines[-1].endswith('\\\\'):
        stat_lines[-1] = stat_lines[-1][:-2]
    latex_stats_without_header_footer = '\n'.join(stat_lines)

    return latex_stats_without_header_footer


def summary_statistics_to_latex(df, outfname):
    # df = pd.read_csv('D:\doing\dc_data_2014.csv')
    sum_stat = df.describe().transpose()
    sum_stat = sum_stat.drop(columns=['25%','50%','75%'])
    sum_stat = sum_stat[sum_stat.index.isin(list(convert_names.keys()))]
    sum_stat = sum_stat.reindex(list(convert_names.keys()))
    sum_stat = sum_stat.rename(index=convert_names)
    sum_stat['unit'] = sum_stat.index.map(units)
    sum_stat = sum_stat.reset_index()
    sum_stat = sum_stat[list(col_names.keys())]
    sum_stat = sum_stat.rename(columns=col_names)
    sum_stat = sum_stat.dropna(how='any', axis=0)

    sums_latex = sum_stat.to_latex(index=False, float_format='%.3f')
    lines = sums_latex.split('\n')[4:-3]
    latex_sumsi = '\n'.join(lines)

    with open(outfname, "w") as tex_file:
        tex_file.write(latex_sumsi)
