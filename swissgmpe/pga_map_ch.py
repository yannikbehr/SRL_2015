#!/usr/bin/env python
"""
Using the GMPE of ?? for Switzerland generate
a map of pga for different magnitudes.
Created on Feb 18, 2014

@author: behry
"""

import os
import sys
sys.path.append('/home/behry/mypy')
import progressbar as pg
import numpy as np
import pyproj
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.pyplot import cm
from mpl_toolkits.basemap import Basemap
from scipy.stats import scoreatpercentile
import json
from swissgmpe import SwissGmpe


def pgamap(mag, latmin, latmax, lonmin, lonmax, ngp, tlat, tlon):
    sg = SwissGmpe()
    lat = np.linspace(latmin, latmax, ngp)
    lon = np.linspace(lonmin, lonmax, ngp)
    sg.grid_setup(lat, lon)
    longrd, latgrd = np.meshgrid(lon, lat)
    # Distinction between alpine and foreland regions
    thresh = 0.39 * longrd + 44

    forelandmask = np.where(latgrd > thresh, 1, 0)
    alpinemask = np.where(latgrd > thresh, 0, 1)
    g = pyproj.Geod(ellps='sphere')
    az, baz, dist = g.inv(longrd, latgrd, np.ones((ngp, ngp)) * tlon,
                          np.ones((ngp, ngp)) * tlat)
    r, pgaforeland = sg.get_pga(mag, dist / 1000., 'foreland', 'median', sn=True)
    r, pgaalpine = sg.get_pga(mag, dist / 1000., 'alpine', 'median', sn=True)
    pga = pgaforeland * forelandmask + pgaalpine * alpinemask
    sg.get_amplification()
    return lon, lat, pga

if __name__ == '__main__':
    tlat = 47.552107
    tlon = 8.228492
    # Number of gridpoints in either direction
    ngp = 201
    mag = 6.25
    latmin, latmax, lonmin, lonmax = 45.4, 48.3, 5.6, 11.1
    lon, lat, pga = pgamap(mag, latmin, latmax, lonmin, lonmax, ngp, tlat, tlon)
    boxout = (45, 48.5, 5, 12, 47)
    meridians = np.arange(5, 12, 2)
    parallels = np.arange(44, 49, 2)
    llcrnrlat, urcrnrlat, llcrnrlon, urcrnrlon, lat_ts = boxout
    m = Basemap(projection='merc', llcrnrlat=llcrnrlat,
                urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                urcrnrlon=urcrnrlon, lat_ts=lat_ts, resolution='i')

    xmin, ymin = m(lonmin, latmin)
    xmax, ymax = m(lonmax, latmax)
    nx = int((xmax - xmin) / 5000.) + 1
    ny = int((ymax - ymin) / 5000.) + 1
    dat, x, y = m.transform_scalar(pga / 1000., lon, lat, nx, ny,
                                   returnxy=True, masked=True)

    # m.fillcontinents(zorder=0)
    cmap = cm.get_cmap('jet')
    cf = m.contourf(x, y, dat, cmap=cmap)
    cs = m.contour(x, y, dat, colors='black',
                   levels=[0.01, 0.02, 0.04, 0.06, 0.1],
                   linestyles='solid', linewidths=3)
    plt.clabel(cs, fmt="%.2f ", fontsize=12, colors='black')

    m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='lightgray',
                    linewidth=0.5, zorder=0)
    m.drawparallels(parallels, labels=[1, 0, 0, 0], color='lightgray',
                    linewidth=0.5, zorder=0)
    m.drawcoastlines(zorder=2)
    m.drawcountries(linewidth=1.0, zorder=2)
    line = 0.32 * np.arange(5, 12, 0.1) + 44.4
    x, y = m(np.arange(5, 12, 0.1), line)
    # m.plot(x, y, 'k-')
    plt.colorbar(cf)
    plt.show()
