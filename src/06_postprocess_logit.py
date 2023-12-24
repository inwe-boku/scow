# %% imports
import os
import logging
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import t
import matplotlib.pyplot as plt
import seaborn as sns
from config import ROOTDIR
from logging_config import setup_logging

setup_logging()

# %% settings
ttype = 'min_lcoe'
sc_name = ['cost']
years = ['2016']  # ['2014', '2016']
mnames = ['base']  # , 'full']
campaigns = [f'{sc}_{yr}_{mn}' for sc in sc_name for yr in years for mn in mnames]

# leaf_type_dummy = True  # set to True if "leaf type" is coded as dummy for "broadleaved" and "coniferous".

# %% retrieve coefficients, stats and implied valuation coefficients
results = {}
#for mname in mnames:
for camp in campaigns:

    ests = pd.read_csv(ROOTDIR / f'data/results/lgtr_coefs_{camp}.csv')
    n = len(ests.groupby('mod').count())

    if ests['term'].str.contains('scalePar').any():
        ests.loc[ests['term'].str.contains('scalePar'), 'term'] = 'lcoe'

    coefs = ests[['term', 'estimate']].groupby(['term']).mean()
    coefs['std'] = ests[['term', 'estimate']].groupby(['term']).std()
    coefs['tstat'] = coefs['estimate'] / coefs['std']
    coefs['pval'] = t.sf(np.abs(coefs['tstat']), n-1)
    coefs['valuation'] = coefs['estimate'] / coefs.loc[f'lcoe_{ttype}', 'estimate']
    coefs.loc[coefs['pval'] <= 0.01, 'significance'] = '***'
    coefs.loc[(coefs['pval'] > 0.01) & (coefs['pval'] <= 0.05), 'significance'] = '**'
    coefs.loc[(coefs['pval'] > 0.05) & (coefs['pval'] <= 0.1), 'significance'] = '*'
    coefs.loc[coefs['pval'] > 0.1, 'significance'] = ''

    loglik = pd.read_csv(ROOTDIR / f'data/results/lgtr_loglik_{camp}.csv')

    results[camp] = {
        'coefs': coefs,
        'loglik': loglik,
    }
    logging.info(f'Model report for {camp} generated.')
    del ests, n

# %% compute marginal social cost of each variable and total social cost
for camp in campaigns:
    sc, yr, _ = camp.split('_')
    geodat = xr.open_dataset(ROOTDIR / f'data/preprocessed/geodat_Niederoesterreich_touch_{sc}_{yr}.nc')

    #for key in results.keys():
    soco = {}
    local_social_cost = geodat['prx_wohnwidmung'] * 0

    for var in list(results[camp]['coefs'].index):
        cost_template = geodat['prx_wohnwidmung'] * 0

        if (':' not in var) and (ttype not in var):
            soco[var] = cost_template + (results[camp]['coefs'].loc[var, 'valuation'] * geodat[var])
        elif (':' not in var) and (ttype in var):
            varshort = var.removesuffix(f'_{ttype}')
            soco[var] = cost_template + (results[camp]['coefs'].loc[var, 'valuation'] * geodat[varshort].sel(turbine_models=ttype))
        else:
            interaction = var.split(':')
            soco[var] = cost_template + (results[camp]['coefs'].loc[var, 'valuation'] * geodat[interaction[0]] * geodat[interaction[1]])

        if 'lcoe' not in var:
            local_social_cost = local_social_cost + soco[var]

    local_social_cost.name = "Local Social Cost"
    local_social_cost.attrs['unit'] = "€/MWh"
    # scale local social cost so that wind turbines are not beneficial locally
    if local_social_cost.min() < 0:
        local_social_cost = local_social_cost - local_social_cost.min()
    soco['Local Social Cost'] = local_social_cost
    social_cost = xr.Dataset(soco)
    social_cost.to_netcdf(ROOTDIR / f'data/results/social_cost_{camp}.nc')
    logging.info(f'Social cost computed and saved to: {ROOTDIR}/data/results/social_cost_{camp}.nc')
    del soco, local_social_cost, var, cost_template, varshort, interaction

# %% plot map of local social cost and LCOE
lcoe = geodat['lcoe']
lcoe.name = "Quasi-LCOE"
lcoe.attrs['unit'] = '€/MWh'
if not os.path.exists(ROOTDIR / 'figures'):
    os.mkdir(ROOTDIR / 'figures')
    logging.info(f'Created directory {ROOTDIR}/figures')

for camp in campaigns:
    social_cost = xr.open_dataset(ROOTDIR / f'data/results/social_cost_{camp}.nc')
    loc_soc_cost = social_cost['Local Social Cost']
    loc_soc_cost.attrs['unit'] = '€/MWh'

    fig, axs = plt.subplots(figsize=(7, 5))
    loc_soc_cost.plot(ax=axs, robust=True)  #vmin=80, vmax=160)  #, vmin=50, vmax=200)
    # lcoe.sel(turbine_models=key).plot(ax=axs[1], robust=True)  # vmin=40, vmax=120)  #, vmin=50, vmax=200)
    axs.set_title("")
    axs.axis('off')
    axs.set_title("")
    axs.axis('off')
    plt.tight_layout()
    plt.savefig(ROOTDIR / f'figures/map_soco_{camp}.png', dpi=200)
    plt.close()

    fig, axs = plt.subplots(figsize=(7, 5))
    # loc_soc_cost.plot(ax=axs[0], robust=True)  #vmin=80, vmax=160)  #, vmin=50, vmax=200)
    lcoe.sel(turbine_models=ttype).plot(ax=axs, robust=True)  # vmin=40, vmax=120)  #, vmin=50, vmax=200)
    axs.set_title("")
    axs.axis('off')
    axs.set_title("")
    axs.axis('off')
    plt.tight_layout()
    plt.savefig(ROOTDIR / f'figures/map_lcoe_{camp}.png', dpi=200)
    plt.close()

logging.info('Maps of Local Social Cost and LCOE plotted')

# %% plot distributions of LCOE, local social cost and total cost
cut_off = 0.999  #0.999

cost_by_turbine = pd.DataFrame()
for camp in campaigns:
    social_cost = xr.open_dataset(ROOTDIR / f'data/results/social_cost_{camp}.nc')

    cost_distribution = pd.DataFrame(data=social_cost['Local Social Cost'].to_pandas().stack().reset_index(drop=True), columns=['Cost [€/MWh]'])
    cost_distribution['target'] = 0
    cost_distribution['Cost Component'] = f'Local Social Cost'

    lcoe_at_all = pd.DataFrame(data=social_cost[f'lcoe_{ttype}'].to_pandas().stack().reset_index(drop=True), columns=['Cost [€/MWh]'])
    lcoe_at_all['target'] = 1
    lcoe_at_all['Cost Component'] = f'Quasi-LCOE'  # f'LCOE_{key}'
    # drop LCOE outliers
    lcoe_at_all = lcoe_at_all.loc[lcoe_at_all['Cost [€/MWh]'] <= lcoe_at_all['Cost [€/MWh]'].quantile(cut_off), :]
    cost_distribution = pd.concat([cost_distribution, lcoe_at_all])

    loco = social_cost['Local Social Cost'] + social_cost[f'lcoe_{ttype}']
    loco = pd.DataFrame(data=loco.to_pandas().stack().reset_index(drop=True), columns=['Cost [€/MWh]'])
    # drop total local social cost outliers
    loco_at_all = loco.loc[loco['Cost [€/MWh]'] <= loco['Cost [€/MWh]'].quantile(cut_off), :]
    # loco_at_all = pd.DataFrame(data=(social_cost['Local Social Cost'] + social_cost[f'lcoe_{key}']).to_pandas().stack().reset_index(drop=True), columns=['Cost'])
    loco_at_all['target'] = 2
    loco_at_all['Cost Component'] = f'Total Local Social Cost'  #f'Total SC_{key}'
    cost_distribution = pd.concat([cost_distribution, loco_at_all])
    cost_by_turbine = pd.concat([cost_by_turbine, cost_distribution])

# cost_by_turbine.loc[cost_by_turbine['Cost'] > 250, 'Cost'] = np.nan
# cost_by_turbine = cost_by_turbine.dropna(axis=0, how='any')

    #sns.displot(data=cost_distribution, x='Cost', hue='Cost Component', fill=True, kind='kde', height=5, aspect=1.5,
    #            palette=sns.color_palette('bright')[:3])

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.kdeplot(data=cost_by_turbine, x='Cost [€/MWh]', hue='Cost Component', fill=True)
    plt.grid(alpha=0.4)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(ROOTDIR / f'figures/dist_soco_{camp}.png', dpi=200)
    plt.close()

logging.info('Cost distributions plotted')

# %% plot social cost components and total local social cost
"""
for camp in campaigns:
    social_cost = xr.open_dataset(ROOTDIR / f'data/results/social_cost_{camp}.nc')
    row, col = 4, 6
    fig, axs = plt.subplots(row, col, figsize=(20, 8))
    n = 0
    for var in list(results[camp]['coefs'].index):
        social_cost[var].plot(ax=axs[int(np.floor(n / col)), n % col], robust=True)
        axs[int(np.floor(n / col)), n % col].set_title(var)
        axs[int(np.floor(n / col)), n % col].axis('off')
        n += 1
    total_social_cost = social_cost['Local Social Cost'].plot(ax=axs[-1, -1], robust=True)
    axs[-1, -1].set_title("Total Local Social Cost")
    axs[-1, -1].axis('off')
#    axs[-1, -2].axis('off')
#    axs[-1, -3].axis('off')
#    axs[-1, -4].axis('off')
    plt.tight_layout()
    plt.savefig(ROOTDIR / f'figures/cost_components_{camp}.png', dpi=200)
    plt.close()

logging.info('Cost components plotted')
"""
