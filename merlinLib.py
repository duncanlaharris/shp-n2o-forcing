# useful functions and variables for MERLIN analysis
import functools
import pathlib

import cf_units
import iris
import iris.analysis
import iris.coord_categorisation
import iris.exceptions
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

clean_cache = True  # set True should you want to clean the cache out
debug = True  # set True for messages
fsize=[11.7, 8.3] # A4 in inches -- curses to matplotlib for using inches..
_merlin_root_dir = pathlib.Path("/users/bedunkadunk/documents/SHP")
#_merlin_root_dir = pathlib.path("c:/users/stett2/data/merlin/FAMOUS")

land_frac_path = _merlin_root_dir/'ancil'/'qrparm.mask_frac.nc'
land_frac= iris.load_cube(str(land_frac_path))
def file_path(experiment, mean_dir, variable):
    """
    Work out path to variable given experiment, variable and mean_dir
    """


    path = _merlin_root_dir / experiment / 'agg' / mean_dir / (experiment + '-' + mean_dir + '-' + variable)
    if not path.exists():
        print("Path ", path, " Does not exist")
    return path


def ocn_vol(file):
    """
    Compute the volume of the ocean from first time in a file
    """
    cube = iris.load_cube(str(file))[0, ...]  # get in file and extract first time.
    cube.data = cube.data / cube.data  # make it all 1 (where not masked)
    layer_thickness = cube.coord('depth').bounds[:, 1] - cube.coord('depth').bounds[:, 0]
    grid_areas = iris.analysis.cartography.area_weights(cube)
    volume_weights = layer_thickness[:, np.newaxis, np.newaxis] * grid_areas
    vol = cube.collapsed(['depth', 'latitude', 'longitude'], iris.analysis.SUM, weights=volume_weights)
    return vol.data


@functools.lru_cache(maxsize=6000)
def read_land(var, exper, use_cache=True):
    if var == 'TotalLandC':
        vars = ['TotalSoilC', 'TotalVegC']
    elif var == 'TotalLandC_ZM':
        vars = ['TotalSoilC_ZM', 'TotalVegC_ZM']
    else:
        raise NotImplementedError(f"Not implemented {var}")

    cube = 0
    for v in vars:
        cube += read_data(v, exper, use_cache=use_cache)
        print(v,cube.data.max())
    return cube
# variables for human editing..
# these might be be better in a separate file...
lookup = pd.read_csv('lookup.csv', header=0, index_col=0)
print(lookup.type)
references = {'Historical': 'xnyvh', 'Control': 'xoauj', 'Spinup': 'xnnqm'}  # reference experiments

_ocn_file = file_path(references['Historical'], 'opx', 'sea_water_potential_temperature.nc')
ocnVol = ocn_vol(str(_ocn_file))  # volume of ocean
mAtmCol = 98630.1 / 9.8  # mass of column (Kg) of atmosphere/m^2 -- int of pstar for one exper/time -- but model should conserve
MoAtm = ((
                     4 * np.pi * 6371229.0 ** 2 * 98630.1) / 9.8) * 1e3  # Mass of atmosphere in grams - 6371229.0 = radius earth in m, 98630.1 = surface pressure in hPa
areaEarth = 4 * np.pi * 6371229.0 ** 2  # area of earth in m^2
# if var_lookup is changed then need to remove the cached data in the data directory (or at least the ones that were affected by the change
PgC = cf_units.Unit('Pg C')  # units for Pg Carbon
PgC_per_yr= cf_units.Unit('Pg/yr C')
## THINK ocean carbon is mol/liter of water = 12.0107 g/liter = 12e-3 g/m^3 = 12e-18 Pg/m^2
ocn_C_scale = 12.0107e-18 # mass of one mole of C in Pg...
CO2_to_C = 12.0107/44.0095 # convert mass CO2 to mass C.
secs_year = 360.*24.*60.*60. # number of seconds in UM year.
var_lookup = {
    'VAOT': dict(file='sea_water_potential_temperature.nc', volAvg=True, Ocean=True),
    'OcnT': dict(file='sea_water_potential_temperature.nc', Ocean=True, thick_scale=True, method=iris.analysis.SUM,
                 scale=1.0 / ocnVol),
    'ocnT_ZM': dict(file='sea_water_potential_temperature.nc', Ocean=True, thick_scale=True, method=iris.analysis.SUM,
                    scale=1.0 / ocnVol, collapse_dims=['longitude']),
    'SST': dict(file='sea_water_potential_temperature.nc', Ocean=True,
                Constraint=iris.Constraint(model_level_number=1)),
    # CO2 emissions -- Kg CO2/m^2/second.
    'CO2emis': dict(file='UM_m01s00i251_vn405.0.nc',Land=True,
                    method=iris.analysis.SUM,scale=CO2_to_C*secs_year*1e-12, units=PgC_per_yr),
    'CO2atm': dict(file='mass_fraction_of_carbon_dioxide_in_air.nc',
                   Constraint=iris.Constraint(model_level_number=11), method=iris.analysis.SUM,
                   scale=CO2_to_C * mAtmCol * 1e-12, units=PgC),
    'MLD': {'file': 'UM_m02s00i137_vn405.0.nc', 'Ocean': True},
    'SAT': {'file': 'surface_temperature.nc'},

    'ocCarbon': dict(file='UM_m02s00i103_vn405.0.nc', Ocean=True, thick_scale=True, scale=ocn_C_scale,
                     method=iris.analysis.SUM, units=PgC),
    'ocCarbonZM': dict(file='UM_m02s00i103_vn405.0.nc', Ocean=True,
                       collapse_dims=['longitude'], thick_scale=True, scale=ocn_C_scale, method=iris.analysis.SUM,
                       units=PgC),
    'TotalOcC': dict(file='UM_m02s00i103_vn405.0.nc', Ocean=True, scale=ocn_C_scale, volAvg=True,
                     method=iris.analysis.SUM, units=PgC),
    'TotalVegC': dict(file='vegetation_carbon_content.nc', Land=True,
                      method=iris.analysis.SUM, scale=1e-12, units=PgC),
    # total veg carbon as Pg. Land Sfc values are kg/m^2
    'TotalSoilC': dict(file='soil_carbon_content.nc', Land=True,
                       method=iris.analysis.SUM, scale=1e-12, units=PgC),
    'TotalLandC':dict(func=read_land),
    # total veg carbon as Pg
    'Litter': dict(file='UM_m01s19i005_vn405.0.nc', Land=True,
                         method=iris.analysis.SUM, scale=1e-12, units=PgC_per_yr),
    'PlantResp':dict(file='plant_respiration_carbon_flux.nc', Land=True,
                         method=iris.analysis.SUM, scale=1e-12*secs_year, units=PgC_per_yr),
    'SoilResp': dict(file='soil_respiration_carbon_flux.nc', Land=True,
                            method=iris.analysis.SUM, scale=1e-12 * secs_year, units=PgC_per_yr),
    # flux as Pg per year.
    'TotalVegC_ZM': dict(file='vegetation_carbon_content.nc', Land=True,
                         method=iris.analysis.SUM, scale=1e-12,
                         units=PgC, collapse_dims=['longitude']),
    # total veg carbon as Pg. Land Sfc values are kg/m^2
    'TotalSoilC_ZM': dict(file='soil_carbon_content.nc', Land=True,
                          method=iris.analysis.SUM, scale=1e-12, units=PgC,
                          collapse_dims=['longitude']),
    # total veg carbon as Pg
    'TotalLitterC_ZM': dict(file='UM_m01s19i005_vn405.0.nc', Land=True,
                            method=iris.analysis.SUM,
                            scale=1e-12*secs_year, units=PgC_per_yr, collapse_dims=['longitude']),
    # total litter flux as Pg per year
    'TotalOcC_ZM': dict(file='UM_m02s00i103_vn405.0.nc', Ocean=True, scale=ocn_C_scale,
                        thick_scale=True, method=iris.analysis.SUM, units=PgC, collapse_dims=['longitude', 'depth']),
    'icef': dict(file='sea_ice_area_fraction.nc', Ocean=True, method=iris.analysis.SUM, scale=1e-12),  # as million km^2
    # ice area in milion km^2
    'NPP': dict(file='net_primary_productivity_of_carbon.nc', Land=True,
                scale=1e-12 * secs_year,
                method=iris.analysis.SUM, units=PgC_per_yr),  # convert to PgC/year
    ## non conservation term. Unit is microMol/liter/day. Thickness is 10m. So sfc layer 10m^3 = 1e3 liters.
    ## ignore -- is actually milli-moles per liter per second (or moles per litre over the whole of the 10m top layer)
    'nonCons':dict(file='UM_m02s30i292_vn405.0.nc',Ocean=True,method=iris.analysis.SUM,
                   scale = ocn_C_scale*secs_year*10, units=PgC_per_yr),

}


@functools.lru_cache(maxsize=6000)
def cache_constraint(**kwargs):
    """
    Cache safe interface to iris.Constraint. Simply passes all arguments in.
    """
    raise NotImplemented
    return iris.Constraint(**kwargs)


def var_properties(var):
    """
    Return properties for specified climate variable. Uses var_lookup
    """

    if var not in var_lookup:
        raise Exception(f"Do not know about {var}")
    if var_lookup[var].get('Ocean', False):
        dir = 'opy'
    else:
        dir = 'apy'

    scale = var_lookup[var].get('scale', 1.0)
    volAvg = var_lookup[var].get('volAvg', False)
    file = var_lookup[var].get('file')
    constraint = var_lookup[var].get('Constraint')
    method = var_lookup[var].get('method', iris.analysis.MEAN)
    thick_scale = var_lookup[var].get('thick_scale', False)
    collapse_dims = var_lookup[var].get('collapse_dims', ['longitude', 'latitude'])  # TODO merge volAvg into this logic
    units = var_lookup[var].get('units')
    if var_lookup[var].get('Land',False):
        lf = np.squeeze(land_frac.data)
    else:
        lf =1.0

    func = var_lookup[var].get('func',None)
    return dict(dir=dir, volAvg=volAvg, file=file, constraint=constraint, scale=scale,
                method=method, thick_scale=thick_scale, lf=lf,func=func,
                collapse_dims=collapse_dims, units=units)


def addAuxCoords(cube):
    """
    Add useful aux coords to a cube
    :param cube: cube to be modified
    :return: nada as cube modified in place
    """
    cube.coord('longitude').circular = True  # make longitude circular
    try:
        iris.coord_categorisation.add_year(cube, 'time')  # add year
        # iris.coord_categorisation.add_month(cube, 'time')  # add month
        iris.coord_categorisation.add_month_number(cube, 'time')  # add month
    except (iris.exceptions.CoordinateNotFoundError, ValueError):
        pass
    for bndCoord in ['time', 'longitude', 'latitude']:
        try:
            cube.coord(bndCoord).guess_bounds()
        except (ValueError, iris.exceptions.CoordinateNotFoundError):
            pass


# ------------
def properties(series,**override):
    ref_lookup = {'Historical': {'colour': 'red', 'linestyle': 'solid'},
                  'Control': {'colour': 'blue', 'linestyle': 'solid'},
                  }
    colour_lookup = dict(Historical={1000: 'firebrick', 500: 'red', 250: 'coral'},
                         Control={1000: 'blue', 500: 'cornflowerblue', 250: 'cyan'})
    alpha_lookup = {1000: 1.0, 500: 0.5, 250: 0.25}
    marker_lookup = {50: 'd', 100: 'h', 200: '*'} # indexed by time

    colour = ref_lookup[series.Reference]['colour']  #
    colour = colour_lookup[series.Reference][1000]
    # alpha = alpha_lookup[series.Carbon]
    alpha = 0.8
    marker = marker_lookup[50]
    result = dict(color=colour, alpha=alpha, marker=marker, linewidth=2)
    result.update(**override)
    return result


def proc_all(var, plot_constraint=None):
    """
    Process all data for a specified variable.

    Args:
        var: str, specidies the variable
        plot_constraint: Iris constraint object, constrains data specified by var
    Returns:
        refs: dict of Iris cubes for keys in global variable references
        diffs: dict of Iris cubes for keys in lookup.csv which access values
            of difference between control and a particular experiment run
        timeseries: dict of Iris cubes for keys in lookup.csv which access
            values of corresponding experiment runs
    """

    refs = {}
    for k, exper in references.items():
        # move the TotalLandC logic into read_data.

        if var is 'Forcing':
            refs[k] = read_data('CO2atm', exper)
        else:
            refs[k] = read_data(var, exper)

    diffs = {}
    timeseries = {}
    for exper, series in lookup.iterrows():

        if var is 'Forcing':  # compute the forcing
            ts = read_data('CO2atm', exper)
            t = diff(ts, refs[series.Reference], ratio=True)
            t = 5.35 * iris.analysis.maths.log(t)  # Myhre formula.
            ts = t

        else:
            ts = read_data(var, exper)
            t = diff(ts, refs[series.Reference])
        diffs[exper] = t.extract(plot_constraint)
        print(f"Processed {exper} for var: {var}")
        timeseries[exper] = ts
    return (refs, diffs, timeseries)


def read_file(var,exper):
    """
    Read file for given variable and experiment
    """
    varprops = var_properties(var)  # properties for variable
    dir = varprops.pop('dir')
    file = varprops.pop('file')
    path = file_path(exper, dir, file)  # path
    cube = iris.load_cube(str(path))  # read the  data
    addAuxCoords(cube)

    return cube

@functools.lru_cache(maxsize=6000)
def read_data(var, exper, use_cache=True):
    """
    Read and process data for variable and experiment
    
    Args:
        var: str
        exper: str
        use_cache: bool
    Returns:
        cube: Iris cube, contains data specified by var in space and time
    """
    # various magic things..


    varprops = var_properties(var)  # properties for variable
    func = varprops.pop('func')
    if func is not None:

        cube = func(var,exper,use_cache=use_cache)
        return cube

    dir = varprops.pop('dir')
    file = varprops.pop('file')
    path = file_path(exper, dir, file)  # path
    proc_path = pathlib.Path('data_cache') / ("_".join((var, exper)) + ".nc")

    if (not clean_cache) and use_cache and proc_path.exists():
        if debug:  print("Using cache at ", proc_path)
        cube = iris.load_cube(str(proc_path))  # read the processed data

    else:

        if debug:  print("Computing using ", path, varprops)
        cube = aggLoad(str(path), **varprops)
        if use_cache:
            iris.save(cube, str(proc_path))  # save the data

    return cube


def aggLoad(filepath, volAvg=False, method=iris.analysis.MEAN, constraint=None, lf=1,
            scale=None, thick_scale=False, collapse_dims=['longitude', 'latitude'], units=None,
            collapse=True):  # default method is global mean

    """
    load and do some processing on data
    """

    cube = iris.load_cube(filepath)
    addAuxCoords(cube)
    cube = cube.extract(constraint)*lf
    weights = iris.analysis.cartography.area_weights(cube)  # assuming long/lat here.
    dims = collapse_dims[:]
    if volAvg or thick_scale:
        layer_thickness = cube.coord('depth').bounds[:, 1] - cube.coord('depth').bounds[:, 0]
        weights = layer_thickness[:, np.newaxis, np.newaxis] * weights
    if volAvg:
        dims.append('depth')

    if collapse:
        mod_global = cube.collapsed(dims, method, weights=weights)
    else:
        mod_global = cube.copy()
        mod_global.data *= weights

    if scale is not None:
        mod_global *= scale

    if units is not None:
        mod_global.units = units

    # mod_global.var_name = var
    return mod_global


def diff(ts, ref, ratio=False):
    """
    Compute difference between ts and reference field.
    """
    interp = ref.interpolate([('time', ts.coord('time').points)],
                             iris.analysis.Linear())
    try:
        if ratio:
            diff = ts / interp
        else:
            diff = (ts - interp)
    except ValueError:  # iris sucks and throws error for this with no useful explanation...
        diff = ts.copy()
        if ratio:
            diff.data = ts.data / interp.data
        else:
            diff.data = ts.data - interp.data

    return diff


def std_error(cube, window=None):
    """
    Compute SD for time-filtered data
    :param cube -- cube (of reference data) to compute sd
    :param window -- window length -- if specified
    """

    c = cube
    if window is not None:
        c = c.rolling_window('time', iris.analysis.MEAN, window)
    # fit 2nd order polynomial.
    ZeroRef_poly, resid_var, rank, singular_values, rcond = np.polyfit(c.coord('time').points, c.data, 2, full=True)
    # resid_var is sum of squares. Need it as an average...
    var = resid_var / c.data.shape[0]  # assuming time is 1st dim
    sd = np.sqrt(var)
    return sd


@functools.lru_cache(maxsize=6000)
def read_cube(var, exper, start_yr, end_yr, weight=False):
    varprops = var_properties(var)  # properties for variable
    dir = varprops.pop('dir')
    file = varprops.pop('file')
    scale = varprops.pop('scale')
    units = varprops.pop('units')
    path = file_path(exper, dir, file)  # path
    cube = iris.load_cube(str(path))
    addAuxCoords(cube)
    cube = cube.extract(iris.Constraint(year=lambda cell: start_yr <= cell < end_yr))
    # cube = cube.collapsed('time', iris.analysis.MEAN)
    if weight:  # weight the variable -- usually to convert to totals.
        weights = iris.analysis.cartography.area_weights(cube)  # assuming long/lat here.
        if varprops['volAvg'] or varprops['thick_scale']:
            layer_thickness = cube.coord('depth').bounds[:, 1] - cube.coord('depth').bounds[:, 0]
            weights = layer_thickness[:, np.newaxis, np.newaxis] * weights
        cube = cube * weights

    if scale is not None:
        cube *= scale

    if units is not None:
        cube.units = units
    cube.var_name = var
    return cube


## make a label object
class plotLabel:
    """
    Class for plotting labels on sub-plots
    """

    def __init__(self, upper=False, roman=False):
        """
        Make instance of plotLabel class
        parameters:
        :param upper -- labels in upper case if True
        :param roman -- labels use roman numbers if True
        """

        import string
        if roman:  # roman numerals
            strings = ['i', 'ii', 'iii', 'iv', 'defaultCov', 'vi', 'vii', 'viii', 'ix', 'x', 'xi', 'xii']
        else:
            strings = [x for x in string.ascii_lowercase]

        if upper:  # upper case if requested
            strings = [x.upper() for x in strings]

        self.strings = strings
        self.num = 0

    def label_str(self):
        """
        Return the next label
        """
        string = self.strings[self.num] + " )"
        self.num += 1
        self.num = self.num % len(self.strings)
        return string

    def plot(self, ax=None, where=None):
        """
        Plot the label on the current axis
        """

        if ax is None:
            plt_axis = plt.gca()
        else:
            plt_axis = ax

        text = self.label_str()
        if where is None:
            x = -0.03
            y = 1.03
        else:
            (x, y) = where

        ax.text(x, y, text, transform=ax.transAxes,
                horizontalalignment='right', verticalalignment='bottom')


def saveFig(fig, name=None, savedir=pathlib.Path("figures"), figtype=None, dpi=None):
    """
    :param fig -- figure to save
    :param name (optional) set to None if undefined
    :param savedir (optional) directory as a pathlib.Path to save figure to . Default is figures path
    :param figtype (optional) type of figure. (If nto specified then png will
    """

    defFigType = '.png'
    if dpi is None:
        dpi = 300
    # set up defaults
    if figtype is None:
        figtype = defFigType
    # work out sub_plot_name.
    if name is None:
        fig_name = fig.get_label()
    else:
        fig_name = name

    outFileName = savedir / (fig_name + figtype)
    fig.savefig(outFileName, dpi=dpi)


def latitudeLabel(value, pos):
    """
    :param values -- value of label
    """

    deg = r'$^\circ$'  # what we need for a degree symbol
    if value < 0:
        end = deg + 'S'
    elif value > 0:
        end = deg + 'N'
    else:
        end = ''

    c = mpl.ticker.ScalarFormatter()
    c.set_scientific(False)
    str = c.format_data_short(abs(value)).strip()  # trailing strip removes whitespace
    str += end
    if abs(value) < 1e-6:
        str = 'Eq'

    return str


lat_format = FuncFormatter(latitudeLabel)
