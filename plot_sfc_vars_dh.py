"""
Reimplementation of vivs code to plot sfc temp, CO2 & Ice.
"""
import iris
import iris.plot
import matplotlib.pyplot as plt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import merlinLib
import numpy as np


plot_constraint = iris.Constraint(year=lambda cell: cell >= 2010)
start_emission = 2020
vars = ['SAT', 'icef', 'NPP','CO2atm','Forcing', 'TotalOcC', 'TotalLandC', 'TotalSoilC']
titles = ['Sfc. Air Temp.', 'Ice Area', 'NPP','Atmospheric CO$_2$', 'Forcing','Ocean Carbon', 'Land Carbon', 'Soil Carbon']
xtitles = ['K', r'$10^{12}$m$^2$', 'Pg a$^{-1}$','Pg C', 'Wm$^{-2}$','Pg C', 'Pg C', 'Pg C']
refs = dict()
diffs = dict()
timeseries = dict()
sd=dict()
for var in vars:
    rr, dd, tt = merlinLib.proc_all(var, plot_constraint=plot_constraint)
    refs[var] = rr
    diffs[var] = dd
    timeseries[var] = tt
    if var is 'Forcing':
        print("SD not well defined for forcing. DO FIX")
        sd[var]=0.0
    else:
        sd[var] = merlinLib.std_error(refs[var]['Control'],window=11)

## plot data
fig, axes = plt.subplots(2, 4, num="sfc", clear=True, figsize=[11.7, 8.3])


for var, ax, ytitle, title in zip(vars, axes.flatten(), xtitles, titles):

    for k, ts in diffs[var].items():
        series = merlinLib.lookup.loc[k]
        prop = merlinLib.properties(series)
        t = ts
        ax.plot(t.coord('year').points, t.data, label=k, markevery=50, **prop)

    ax.set_title(title)
    ax.set_ylabel(ytitle)
    ax.axhline(0.0, linestyle='dashed', color='black')
    ax.fill_between(ax.get_xlim(),2*sd[var],-2*sd[var],color='grey',alpha=0.4)
    # Plot grey bars to show start and end of emissions.
    # assuming annual data..
    y = ax.get_ylim()
    ax.fill_betweenx(y, start_emission, start_emission + 10, alpha=0.6, color='grey')
    for t in merlinLib.lookup.Time.unique():
        ax.fill_betweenx(y, t + start_emission - 10, t + start_emission, alpha=0.4, color='grey')

# done do final figure stuff.
label = merlinLib.plotLabel()  # new set of labels
for ax in axes.flatten():
    label.plot(ax)
fig.tight_layout()
fig.show()
merlinLib.saveFig(fig)

## plot the temp
import functools


@functools.lru_cache(maxsize=6000)
def read_cube(var, exper):
    varprops = merlinLib.var_properties(var)  # properties for variable
    path = merlinLib.file_path(exper, varprops['dir'], varprops['file'])  # path
    cube = iris.load_cube(str(path))
    merlinLib.addAuxCoords(cube)
    cube = cube.extract(iris.Constraint(year=lambda cell: 2450 <= cell < 2500))
    cube = cube.collapsed('time', iris.analysis.MEAN)
    return cube

@functools.lru_cache(maxsize=6000)
def mean_variance(var, exper,ratio=False,crit=None):
    varprops = merlinLib.var_properties(var)  # properties for variable
    path = merlinLib.file_path(exper, varprops['dir'], varprops['file'])  # path
    cube = iris.load_cube(str(path))
    merlinLib.addAuxCoords(cube)
    if ratio:# varn rel to 2450-2500.
        ref_cube=cube.extract(iris.Constraint(year=lambda cell: 2450 <= cell < 2500))
        ref_cube = ref_cube.collapsed('time', iris.analysis.MEAN)

        if crit is not None:
            L=np.abs(ref_cube.data) < crit
            ref_cube.data.mask[L]=True
        cube=cube/ref_cube
    cube=cube.extract(iris.Constraint(year=lambda cell: cell< 2450))
    cube = cube.collapsed('longitude',iris.analysis.MEAN)
    cube = ts
    variance = cube.collapsed('time',iris.analysis.VARIANCE)
    if not ratio:
        variance *= 2 # double variance as difference. Ratio calc includes that.
    return variance

@functools.lru_cache(maxsize=6000)
def mean_dif(var, exper,ratio=False,crit=None):
    cube = read_cube(var, exper)
    ref_exper = merlinLib.references[merlinLib.lookup.loc[exper, 'Reference']]
    cube2 = read_cube(var, ref_exper)
    if ratio:

        delta=cube/cube2
        if crit is not None:
            L=np.abs(cube2.data) < crit
            delta.data.mask[L]=True

    else:
        delta = cube - cube2

    return delta


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

"""
import cartopy.crs as ccrs
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter

diffs=dict()
variance=dict()
for var in ['SAT','NPP']:
    d = {}
    variance[var] = mean_variance(var,merlinLib.references['Control'],ratio=(var == 'NPP'),crit=5e-9)
    sub = merlinLib.lookup.query('Time ==200')

    for exper, series in sub.iterrows():
        d[exper] = mean_dif(var, exper,ratio=(var == 'NPP'),crit=5e-9)
        if var == 'NPP':
            d[exper]=d[exper]-1.0

    diffs[var] = d

## plot now
import numpy as np
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()

fig = plt.figure(num='delta', clear=True, figsize=[11.7, 8.3])
axMn_T = fig.add_subplot(221, projection=ccrs.PlateCarree())
axZM_T = fig.add_subplot(222,sharey=axMn_T)
axMn_NPP = fig.add_subplot(223, projection=ccrs.PlateCarree())
axZM_NPP = fig.add_subplot(224,sharey=axMn_NPP)
levels_t = [-3, -2, -1, -0.5, -0.25, 0.25, 0.5, 1, 2, 3]
levels_NPP = np.array([0.8,0.9,0.95,1.05,1.1,1.2])-1.
for axMn,axZM,var,levels,title in zip([axMn_T,axMn_NPP],[axZM_T,axZM_NPP],['SAT','NPP'],[levels_t,levels_NPP],
                                      ['Temperature','NPP']):
    diff=diffs[var]['xocyc']
    cf=iris.plot.contourf(diff, levels=levels, axes=axMn)
    CS=iris.plot.contour(diff, levels=levels, axes=axMn)
    axMn.clabel(CS,colors='black',fmt='%3.2f')
    axMn.set_xticks([-180,-120,-60,0, 60, 120,180], crs=ccrs.PlateCarree()) # from cartopy example
    # but note the co-ords need to be given in the right order...
    axMn.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    axMn.coastlines()
    axMn.set_title(rf'{title} $\Delta$ (2450-2500) ')
    axMn.set_xlabel('Longitude')
    axMn.set_ylabel('Latitude')
    axMn.xaxis.set_major_formatter(lon_formatter)
    axMn.yaxis.set_major_formatter(lat_formatter)
    cbar = fig.colorbar(cf,  ax=axMn,ticks=levels, orientation='horizontal',
                        extend='both',drawedges=True)

    sd_zm = variance[var]**0.5
    for k, diff in diffs[var].items():
        series = merlinLib.lookup.loc[k, :]
        prop = merlinLib.properties(series)
        c = diff.collapsed('longitude', iris.analysis.MEAN)# * 1000 / series.Carbon
        y = c.coord('latitude').points
        axZM.plot(c.data, y, markevery=2, ms=10, **prop)
    axZM.fill_betweenx(y, (2*sd_zm).data, (-2*sd_zm).data, color='grey', alpha=0.3)
    axZM.axvline(0.0, linestyle='dashed', color='black')
    axZM.set_xlabel(f"{title} ")
    #axZM.set_ylabel('Latitude')
    axZM.set_title(f"Zonal Mean {title} difference")



label = merlinLib.plotLabel()
for ax in (axMn_T, axZM_T,axMn_NPP,axZM_NPP):
    label.plot(ax)


fig.tight_layout()
# adjust posn
for axMn,axZM,var,levels in zip([axMn_T,axMn_NPP],[axZM_T,axZM_NPP],['SAT','NPP'],[levels_t,levels_NPP]):
    mnbounds= axMn.get_position() # l,b,w,h
    zmbounds = axZM.get_position()
# the /2 factor needed to make this work. My guess as to why is that the axis correctly proportional works
# for this plot gives almost a factor of two white space in the plot. Do not use tight_layout -- it breaks this hack.
# THere must be a better way of doing this....
    newBnds = mnbounds.from_bounds(zmbounds.bounds[0],mnbounds.bounds[1],zmbounds.bounds[2],mnbounds.bounds[3])
    axZM.set_position(newBnds)

fig.show()
merlinLib.saveFig(fig)

##
"""