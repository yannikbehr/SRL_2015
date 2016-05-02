#!/usr/bin/env python
"""
Make a plot of P-wave travel time to 6 stations in the network.

Created on Jun 24, 2013

@author: behry
"""
import sys
sys.path.append('./swissgmpe')
from swissgmpe import SwissGmpe
from pga_map_ch import pgamap
import numpy as np
import matplotlib
try:
    matplotlib.use('WXAgg')
    import wx
except:
    print "WX package for Python not installed"
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.pyplot import cm
from mpl_toolkits.basemap import Basemap
from scipy.stats import scoreatpercentile
from copy import copy
import warnings


class AlertTimeMap:

    def __init__(self, resultsfn, mapbox=(45, 48.5, 5, 12, 47),
                 cmapname='RdBu_r'):
        llcrnrlat, urcrnrlat, llcrnrlon, urcrnrlon, lat_ts = mapbox
        self.m = Basemap(projection='merc', llcrnrlat=llcrnrlat,
                         urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon,
                         urcrnrlon=urcrnrlon, lat_ts=lat_ts,
                         resolution='i')
        self.cmapname = cmapname
        self.extend = 'max'
        self.load_data(resultsfn)

    def load_data(self, resultsfn):
        a = np.load(resultsfn)
        self.tstarget = a['tstarget']
        self.ttP = a['ttP']
        self.dims = self.ttP.shape
        if len(self.dims) > 2:
            self.lat = a['lat'][:, :, 0]
            self.lon = a['lon'][:, :, 0]
        else:
            self.lat = a['lat']
            self.lon = a['lon']
        self.x, self.y = self.m(self.lon, self.lat)

    def plot1D(self, fig, ax, scale=True, meridians=np.arange(5, 12, 2),
               parallels=np.arange(44, 49, 2), vmin=0., vmax=15):
        """
        Show alert times for selected events.
        """
        ttPmed = np.median(self.ttP, axis=0)
        cmap = cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax),
                                 cmap=self.cmapname)
        for lt, ln, delay in zip(self.lat, self.lon, ttPmed):
            cl = cmap.to_rgba(delay)
            x, y = self.m(ln, lt)
            self.m.plot(x, y, ms=8, c=cl, marker='o', picker=5.)
        self.m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='lightgray',
                             linewidth=0.5, zorder=0)
        self.m.drawparallels(parallels, labels=[1, 0, 0, 0], color='lightgray',
                             linewidth=0.5, zorder=0)
        self.m.drawcoastlines(zorder=2)
        self.m.drawcountries(linewidth=1.0, zorder=2)
        self.m.drawstates(zorder=2)


    def plot2D(self, fig, ax, depth=8, vs=3.5, nnst=6, procdelay=False,
               scale=False, boxin=(45.4, 48.3, 5.6, 11.1),
               geofilter=None, boxout=(45, 48.5, 5, 12, 47), target=None,
               meridians=np.arange(5, 12, 2), interactive=False,
               parallels=np.arange(44, 49, 2), blindzone=False,
               pga=False, clevels=None, cbshrink=0.8):

            self.m.ax = ax
            if target is not None:
                tln, tlt = target

            latmin, latmax, lonmin, lonmax = boxin
            cmap = cm.get_cmap(self.cmapname)
            cmap.set_over('grey')
            if procdelay:
                self.vmin = 6.
                self.vmax = 25.
                self.cb_label = 'Time since origin time [s]'
                unit = 's'
                clevels_colours = ['lightgray', 'gray']
                dlevel = 0.5
            if target is not None:
                self.extend = 'both'
                self.vmin = -10.
                self.vmax = 60.
                cmap = cm.get_cmap(self.cmapname + '_r')
                self.cb_label = 'Lead time [s]'
                clevels_colours = ['lightgray', 'gray']
                unit = 's'
                dlevel = 0.5
            if not procdelay and target is None:
                self.vmin = 0.
                self.vmax = 15.
                self.cb_label = 'P-wave travel time to %d stations [s]' % nnst
                clevels_colours = ['gray', 'gray']
                unit = 's'
                dlevel = 0.5
            if blindzone:
                self.vmin = 22.
                self.vmax = 55.
                self.cb_label = 'Blind zone [km]'
                clevels_colours = ['lightgray', 'gray']
                unit = 'km'
                dlevel = 0.5
                if pga:
                    self.extend = 'both'
                    self.vmin = 0.01
                    self.vmax = 0.1
                    cmap = cm.get_cmap('PuBu')
                    self.cb_label = 'Maximum pga [g]'
                    unit = 'g'
                    dlevel = 0.005

            # Mask points outside of polygon
            if geofilter is not None:
                if self.dims > 2:
                    rowidx = 0
                    colidx = 0
                    idx = 0
                    ydim = self.lat.shape[1]
                    for _lat, _lon in zip(self.lat.ravel(), self.lon.ravel()):
                        rowidx = idx / ydim
                        colidx = idx - rowidx * ydim
                        idx += 1
                        if not geofilter.point_in_polygon([_lon], [_lat])[0]:
                            self.ttP[:, rowidx, colidx, :] = np.nan
                else:
                    idx = 0
                    for _lat, _lon in zip(self.lat, self.lon):
                        if not geofilter.point_in_polygon([_lon], [_lat])[0]:
                            self.ttP[:, idx] = np.nan
                        idx += 1

            ttPmed = np.median(self.ttP, axis=0)
            ttPmed = np.median(ttPmed, axis=-1)
            print "The minimum alert time is: ", np.ma.masked_invalid(ttPmed).min()
            if target is not None:
                print self.tstarget.min(), self.tstarget.max()
                ttPmed = self.tstarget[:, :, 0] - ttPmed

            if blindzone:
                bz = np.sqrt(ttPmed * ttPmed * vs * vs - depth * depth)
                if pga:
                    sg = SwissGmpe()
                    sg.grid_setup(self.lat, self.lon)
                    mag = 6.5
                    r, pgaforeland = sg.get_pga(mag, bz, 'foreland', 'median')
                    r, pgaalpine = sg.get_pga(mag, bz, 'alpine', 'median')
                    ttPmed = pgaforeland * sg.forelandmask + pgaalpine * sg.alpinemask
                    ttPmed *= 1.7
                    ttPmed /= 1000
                    print ttPmed.min(), ttPmed.max()
                else:
                    ttPmed = bz
                print "The minimum blindzone is: ", ttPmed.min()

            cf = self.m.contourf(self.x, self.y, ttPmed, cmap=cmap,
                                 levels=np.arange(self.vmin,
                                                  self.vmax + dlevel, dlevel),
                                 norm=Normalize(vmin=self.vmin,
                                                vmax=self.vmax),
                                 extend=self.extend)

            if target is not None:
                xt, yt = self.m(tln, tlt)
                self.m.plot(xt, yt, marker='*', ms=14, color='black')

            # Add contour lines
            if clevels is not None:
                for _lev, _col in zip(clevels, clevels_colours):
                    cs = self.m.contour(self.x, self.y, ttPmed,
                                        colors=_col, levels=[_lev],
                                        linestyles='solid', linewidths=3)
                    with warnings.catch_warnings(record=True):
                        plt.clabel(cs, fmt="%d " + unit, fontsize=12, colors=_col)

            if scale:
                cb = fig.colorbar(cf, ax=ax, extend=self.extend,
                                  orientation='vertical',
                                  spacing='uniform', shrink=cbshrink)
                cb.set_label(self.cb_label)
            self.m.drawmeridians(meridians, labels=[0, 0, 0, 1], color='lightgray',
                            linewidth=0.5, zorder=0)
            self.m.drawparallels(parallels, labels=[1, 0, 0, 0], color='lightgray',
                            linewidth=0.5, zorder=0)
            self.m.drawcoastlines(zorder=2)
            self.m.drawcountries(linewidth=1.0, zorder=2)
            self.m.drawstates(zorder=2)

            if pga and target is not None:
                mag = 6.25
                lon, lat, pga = pgamap(mag, latmin, latmax,
                                       lonmin, lonmax, self.ngp, tlt, tln)
                cs = self.m.contour(self.x, self.y, pga / 1000., colors='black',
                               levels=[0.01, 0.02, 0.04, 0.06, 0.1],
                               linestyles='solid', linewidths=3)
                plt.clabel(cs, fmt="%.2f ", fontsize=12, colors='black')
                txtx, txty = self.m(5.5, 45.3)
                ax.text(txtx, txty, 'Magnitude=%.2f + 1$\sigma$' % mag,
                        horizontalalignment='left',
                        verticalalignment='center',
                        backgroundcolor='white')
            return self.m

    def plot_stations(self, fig, ax, networkinfo, interactive=False,
                      networklabels=False):
        # Plot station locations
        dataX = []
        dataY = []
        values = []
        networks = networkinfo.get_networks()
        for net in networks.keys():
            stlat = networks[net]['lat']
            stlon = networks[net]['lon']
            names = networks[net]['nm']
            color = networks[net]['color']
            label = networks[net]['label']
            nwcode = np.array(networks[net]['nw'])
            for _ln, _lt, nm, nt in zip(stlon, stlat, names, nwcode):
                x, y = self.m(_ln, _lt)
                values.append('%s.%s' % (nt, nm))
                dataX.append(x)
                dataY.append(y)
                self.m.plot(x, y, color=color, marker='^', mec='black', ms=6,
                            picker=5.)
            if networklabels:
                self.m.plot(x, y, color=color, marker='^', mec='black', ms=5,
                       label=label, ls='None')
        if networklabels:
                ax.legend(numpoints=1, borderpad=0.2, fontsize='small',
                          markerscale=1.5)

        # add an interactive picker
        if interactive:
            try:
                # add a pop-up window showing the station and its value
                tooltip = wx.ToolTip(tip='')
                tooltip.Enable(False)
                tooltip.SetDelay(0)
                fig.canvas.SetToolTip(tooltip)

                def onMotion(event):
                    line2d = event.artist
                    x = line2d.get_xdata()[0]
                    y = line2d.get_ydata()[0]
                    found = False
                    for i in xrange(len(dataX)):
                        radius = 5
                        if abs(x - dataX[i]) < radius and abs(y - dataY[i]) < radius:
                            tip = '%s' % values[i]
                            tooltip.SetTip(tip)
                            tooltip.Enable(True)
                            found = True
                            break
                    if not found:
                        tooltip.Enable(False)
                fig.canvas.mpl_connect('pick_event', onMotion)
            except Exception, e:
                print "Cannot add wx.ToolTip: ", e

    def error_plot(self, ax_lb, ax_ub, cax, cborientation='vertical'):
        # plot the error map
        ttP_lb = np.zeros((self.dims[1::]))
        ttP_ub = ttP_lb.copy()
        for _i1 in xrange(self.dims[1]):
            for _i2 in xrange(self.dims[2]):
                for _i3 in xrange(self.dims[3]):
                    ttP_lb[_i1, _i2, _i3] = scoreatpercentile(self.ttP[:, _i1, _i2, _i3], 16)
                    ttP_ub[_i1, _i2, _i3] = scoreatpercentile(self.ttP[:, _i1, _i2, _i3], 84)

        mlb = copy(self.m)
        mlb.ax = ax_lb
        mub = copy(self.m)
        mub.ax = ax_ub

        cmap = cm.get_cmap(self.cmapname)
        cmap.set_over('grey')

        mlb.contourf(self.x, self.y, ttP_lb[:, :, 0], cmap=cmap,
                     levels=np.arange(self.vmin, self.vmax + 0.5, 0.5),
                     norm=Normalize(vmin=self.vmin, vmax=self.vmax),
                     extend=self.extend)
        mub.contourf(self.x, self.y, ttP_ub[:, :, 0], cmap=cmap,
                     levels=np.arange(self.vmin, self.vmax + 0.5, 0.5),
                     norm=Normalize(vmin=self.vmin, vmax=self.vmax),
                     extend=self.extend)
        mlb.drawcoastlines(zorder=2)
        mlb.drawcountries(linewidth=1.0, zorder=2)
        mub.drawcoastlines(zorder=2)
        mub.drawcountries(linewidth=1.0, zorder=2)
        cb = ColorbarBase(cax, cmap=cmap, norm=Normalize(vmin=self.vmin,
                                                         vmax=self.vmax),
                         orientation=cborientation, extend=self.extend)
        cb.set_label(self.cb_label)
        return mlb, mub


if __name__ == '__main__':
    pass
