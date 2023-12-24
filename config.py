# %% imports
from pathlib import Path

# %% global settings
ROOTDIR = Path('c:/git_repos/scow')
GAMSDIR = Path('c:/myprogs/GAMS/45')
RDIR = Path('c:/myprogs/R/R-4.2.2')

gams_conf = {
    'gams_model': ROOTDIR / 'opt/location_selection.gms',
    'gdx_input': ROOTDIR / 'opt/input_data.gdx',
    'gdx_output': ROOTDIR / 'opt',
    'gams_exe': GAMSDIR
}

country = 'AUT'

# optimization settings
nslices = 9
nturbines = 'auto'
spacing = 3
max_turb = 8000
min_turb = 1000

# settings for capacity factor simulation and LCOE computation
BIAS_CORRECTION = 0.71197  # 0.85
TURBINE_COST_SHARE = 0.7  # use only turbine cost in quasi-lcoe computation; share: see nrel-report on wind energy 2023
OM_FIXED = 20  # 20  # EUR/kW
OM_VARIABLE = 8  # 26.4  # 0.008 * 1000  # EUR/kWh
DISCOUNT_RATE = 0.05  # 0.05  # 0.03
LIFETIME = 25  # 25
TURBINE_YEARS = [2014, 2016]

# turbines to consider
turbines = {
    #'Dewind.D6.1250': [1250, 92, 64, 2003],  # 'Dewind.D6.1250'
    #'Dewind.D8.2000': [2000, 100, 80, 2003],  # NA
    'Enercon.E40.500': [500, 50, 40, 1998],  # * 'Enercon.E40 5.40.500'
    'Enercon.E40.600': [600, 65, 40, 2002],  # 'Enercon.E40 6.44.600'
    'Enercon.E66.1800': [1800, 85, 66, 2003],  # 'Enercon.E66 18.70.1800'
    'Enercon.E66.2000': [2000, 98, 66, 2004],  # 'Enercon.E66 20.70.2000'
    #'Enercon.E70.1800': [1800, 86, 70, 2006],
    'Enercon.E70.2000': [2000, 86, 70, 2006],  # 'Enercon.E70 E4.2000'
    'Enercon.E70.2300': [2300, 113, 70, 2010],  # 2012
    'Enercon.E82.2300': [2300, 108, 82, 2012],  # 2012 'Enercon.E82 E2.2300'
    'Enercon.E82.3000': [3000, 78, 82, 2014],  # * 2016
    'Enercon.E92.2350': [2350, 104, 92, 2013],  # * 2015 gelistet seit 2012;
    'Enercon.E101.3050': [3050, 135, 101, 2012],  # * 2014 gelistet seit 2012; https://www.wind-turbine-models.com/turbines/130-enercon-e-101
    'Enercon.E115.3000': [3000, 135, 115, 2014],  # * 2016 gelistet seit 11/2013; https://www.wind-turbine-models.com/turbines/832-enercon-e-115-3.000
    'Enercon.E126.3500': [3500, 135, 126, 2013],  # 2016
    'Enercon.E126.7500': [7500, 138, 126, 2013],  # 2012
    'Enercon.E138.3500': [3500, 131, 138, 2016],  # 2016
    'Enercon.E160.5560': [5560, 166, 160, 2016],  # 2016
    'GE.1.5sl': [1500, 85, 77, 2004],  # 'GE.1.5sl.1500'
    # #'NEG.Micon.750.750': 70,
    # #'NEG.Micon.1500.1500': 60,
    'Nordex.N29.250': [250, 50, 30, 1996],
    # # 'Repower.3XM.3200': 128, REpower.3.4M
    #'Repower.M114.3000': [3000, 143, 114, 2013],
    'REpower.MM82.2000': [2000, 100, 82, 2005],
    'REpower.MM92.2000': [2000, 100, 92, 2012],
    # 'Senvion/REpower.S114.3200': 143,  # 'Senvion.3.2M114.3170'
    # #'Siemens.Bonus.1300': 60,
    'Vestas.V44.600': [600, 63, 44, 1997],
    'Vestas.V47.660': [660, 65, 47, 2000],
    'Vestas.V80.2000': [2000, 100, 80, 2003],  #
    'Vestas.V90.2000': [2000, 105, 90, 2007],
    'Vestas.V100.1800': [1800, 100, 100, 2013],  # 2012
    'Vestas.V100.2000': [2000, 105, 100, 2014],  # * 2016 gelistet seit 10/2014
    'Vestas.V112.3000': [3000, 120, 112, 2015],
    'Vestas.V112.3075': [3075, 120, 112, 2015],
    # 'Vestas.V112 a.3000': 140,
    'Vestas.V112.3300': [3300, 94, 112, 2014],  # 2016 gelistet seit 9/2013
    'Vestas.V112.3450': [3450, 84, 112, 2014],  # 2016
    'Vestas.V117.3300': [3300, 130, 117, 2015],  # 2016
    'Vestas.V126.3300': [3300, 137, 126, 2015],  # * 2016 gelistet seit 9/2013
    'Vestas.V126.3450': [3450, 117, 126, 2015],  #
    'Vestas.V136.4200': [4200, 149, 136, 2016],
    'Vestas.V150.4200': [4200, 155, 150, 2016],
    'Vestas.V162.5600': [5600, 166, 162, 2016],
    'Vestas.V162.7200': [7200, 169, 162, 2016],
}
