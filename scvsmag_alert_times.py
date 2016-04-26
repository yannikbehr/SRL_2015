#!/usr/bin/env python
"""
Plot times between origin time and the first alert issued by scvsmag.
Created on Jun 28, 2013

@author: behry
"""
import glob
import os
import math
from obspy import UTCDateTime
from scvsmaglog_report_parser import parser as scvs_parser
import matplotlib
try:
    matplotlib.use('WXAgg')
    import wx
except:
    print "WX package for Python not installed"
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from mpl_toolkits.basemap import Basemap
from scipy.io.netcdf import NetCDFFile as Dataset
from matplotlib.colors import LightSource
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.pyplot import cm
from matplotlib.patches import Wedge
from scipy.stats import scoreatpercentile
from scipy import spatial
import json


class AlertTimes:
    def __init__(self, datdir, datdir2):
        self.del_coord = {}
        self.datdir = datdir
        self.datdir2 = datdir2
        self.delays = []
        self.del_coord_rt = {}

    def gauss(self, x, y):
        fitfunc = lambda p, x: (1 / np.sqrt(2 * np.pi * p[0] ** 2)) * np.exp(-(x - p[1]) ** 2 / (2 * p[0] ** 2))
        errfunc = lambda p, x, y: fitfunc(p, x) - y
        gaussian = lambda m, s, x: (1 / np.sqrt(2 * np.pi * s ** 2)) * np.exp(-(x - m) ** 2 / (2 * s ** 2))
        p0 = [2., 15.]
        p1, success = optimize.leastsq(errfunc, p0[:], args=(x, y))
        sigm, mean = p1
        return round(abs(sigm), 1), round(mean, 1), gaussian(mean, sigm, x)

    def get_delays(self, new=True, delayfn='./data/alert_times.txt'):
        if new:
            for dirn in glob.glob(os.path.join(self.datdir, '*'))[:]:
                evid = os.path.basename(dirn)
                reps = glob.glob(os.path.join(dirn, '*_report.txt'))
                dtls = os.path.join(dirn, evid + '_details.txt')
                if not os.path.isfile(dtls) or len(reps) < 1:
                    print "problems with event %s" % evid
                    continue
                else:
                    f = open(dtls, 'r')
                    f.readline()
                    line = f.readline()
                    a = line.split()
                    true_ot = UTCDateTime(a[4])
                    lat = float(a[0])
                    lon = float(a[1])
                    true_loc = (lat, lon)
                    true_mag = float(a[3])
                    f.close()
                if len(reps) > 1:
                    # print reps
                    min_tdiff = 1e39
                    for _rr in reps:
                        _t = scvs_parser(_rr, trueot=true_ot, trueloc=true_loc,
                                         format='new')
                        mag, lt, ln, dep, ct, ot, td, ddiff, lkh, nstorig, nstmag = _t[0]
                        _tdiff = abs(UTCDateTime(ct) - true_ot)
                        if _tdiff < min_tdiff:
                            min_tdiff = _tdiff
                            cts = _t
                else:
                    cts = scvs_parser(reps[0], trueot=true_ot, trueloc=true_loc,
                                      format='new')
                mags, lats, lons, deps, ct_s, ots, tds, ddiffs, lkh, nos, nms = cts[0]
                if nos > 6:
                    print "Number of stations used for 1st origin: %d" % nos
                    print evid
                    continue
                mage, late, lone, depe, ct_e, ote, tde, ddiffe, lkh, noe, nme = cts[-1]
                del1 = UTCDateTime(ct_s)
                # The playback runs with station delays have to be made
                # in realtime mode. Therefore origin times are not the original
                # ones anymore. To get an idea of the origin time error I compare
                # the origin time error for the same events from historic playbacks.
                reps2 = glob.glob(os.path.join(self.datdir2, evid, '*_report.txt'))
                if len(reps2) > 1:
                    # print reps
                    min_tdiff = 1e39
                    for _rr in reps2:
                        _t = scvs_parser(_rr, trueot=true_ot, trueloc=true_loc,
                                         format='new')
                        mag, lt, ln, dep, ct, ot, td, ddiff, lkh, no, nm = _t[0]
                        _tdiff = abs(UTCDateTime(ct) - true_ot)
                        if _tdiff < min_tdiff:
                            min_tdiff = _tdiff
                            cts2 = _t
                else:
                    cts2 = scvs_parser(reps2[0], trueot=true_ot, trueloc=true_loc,
                                      format='new')
                ot_del = UTCDateTime(cts2[-1, 5]) - UTCDateTime(true_ot)
                print dirn, del1, ot_del, ddiffe
                self.delays.append(del1 - UTCDateTime(ote) + ot_del)
                self.del_coord[evid] = (lon, lat,
                                        del1 - UTCDateTime(ote) + ot_del,
                                        mags - true_mag)
            fh = open(delayfn, 'w')
            json.dump(self.del_coord, fh)
            fh.close()
        else:
            fh = open(delayfn)
            self.del_coord = json.load(fh)
            fh.close()

    def get_delays_rt(self, fn, repdir, delayfn='./data/alert_times_rt.txt',
                      cutoffmag=2.5, new=True):
        if new:
            fh = open(fn)
            for _l in fh.readlines():
                if _l.startswith('#'):
                    continue
                a = _l.split()
                true_lat = float(a[1])
                true_lon = float(a[2])
                true_loc = (true_lat, true_lon)
                true_dep = float(a[3])
                true_mag = float(a[4])
                true_ot = UTCDateTime(a[5])
                bn = a[-2]
                evid = a[-1].rstrip()
                if bn != 'None':
                    repfn = os.path.join(repdir, bn)
                    # only read new format reports that contain the number of
                    # stations used for the origins
                    try:
                        ct = scvs_parser(repfn, trueot=true_ot, trueloc=true_loc,
                                         format='new')
                    except:
                        print repfn
                        continue
                    mag, lt, ln, dep, ct, ot, td, ddiff, lkh, no, nm = ct[0]
                    if no > 6 or mag < cutoffmag:
                        continue
                    else:
                        del1 = UTCDateTime(ct)
                        self.del_coord_rt[evid] = (true_lon, true_lat,
                                                   del1 - true_ot,
                                                   mag - true_mag)
            fh = open(delayfn, 'w')
            json.dump(self.del_coord_rt, fh)
            fh.close()
        else:
            fh = open(delayfn)
            self.del_coord_rt = json.load(fh)
            fh.close()

    def plot_alert_times(self, fout):
        fig = plt.figure(figsize=(9.6, 3))
        ax = fig.add_subplot(111)
        n, bins, patches = ax.hist(self.delays,
                                   bins=np.arange(0, 40, 1.0),
                                   color='blue',
                                   label='All stations', rwidth=1.0,
                                   histtype='bar', normed=False)
        if False:
            x = bins[0:-1] + np.diff(bins)
            y = n
            sigm, mean, gy = self.gauss(x, y)
        med = np.median(self.delays)
        percentile25 = scoreatpercentile(self.delays, 25)
        percentile75 = scoreatpercentile(self.delays, 75)
        err = ((med - percentile25) + (percentile75 - med)) / 2.
        ax.text(0.6, 0.7, 'Median: %.1f s' % (med), horizontalalignment='left',
                transform=ax.transAxes, color='blue')
        ax.text(0.6, 0.63, '25th percentile: %.1f s' % (percentile25), horizontalalignment='left',
                transform=ax.transAxes, color='blue')
        ax.text(0.6, 0.56, '75th percentile: %.1f s' % (percentile75), horizontalalignment='left',
                transform=ax.transAxes, color='blue')
        ax.set_xlabel('Time since origin time [s]')
        ax.set_ylim(0, 16)

        plt.savefig(fout, dpi=300, bbox_inches='tight')

    def background_map(self, ax):
        m = Basemap(projection='merc', llcrnrlat=44.5,
                    urcrnrlat=48.5, llcrnrlon=5, urcrnrlon=12, lat_ts=47,
                    resolution='i', ax=ax)
        m.drawmapboundary(fill_color='lightblue', zorder=0)
        m.fillcontinents(zorder=0)
        etopofn = '/home/behry/uni/data/etopo1_central_europe_gmt.grd'
        etopodata = Dataset(etopofn, 'r')
        z = etopodata.variables['z'][:]
        x_range = etopodata.variables['x_range'][:]
        y_range = etopodata.variables['y_range'][:]
        spc = etopodata.variables['spacing'][:]
        lats = np.arange(y_range[0], y_range[1], spc[1])
        lons = np.arange(x_range[0], x_range[1], spc[0])
        topoin = z.reshape(lats.size, lons.size, order='C')
        # transform to nx x ny regularly spaced 5km native projection grid
        nx = int((m.xmax - m.xmin) / 5000.) + 1; ny = int((m.ymax - m.ymin) / 5000.) + 1
        topodat, x, y = m.transform_scalar(np.flipud(topoin), lons, lats, nx, ny, returnxy=True)
        ls = LightSource(azdeg=300, altdeg=15, hsv_min_sat=0.2, hsv_max_sat=0.3,
                         hsv_min_val=0.2, hsv_max_val=0.3)
        # shade data, creating an rgb array.
        rgb = ls.shade(np.ma.masked_less(topodat / 1000.0, 0.0), cm.gist_gray_r)
        m.imshow(rgb)
        m.drawmeridians(np.arange(6, 12, 2), labels=[0, 0, 0, 1], color='white',
                        linewidth=0.5, zorder=0)
        m.drawparallels(np.arange(44, 50, 2), labels=[1, 0, 0, 0], color='white',
                        linewidth=0.5, zorder=0)
        m.drawcoastlines(zorder=1)
        m.drawcountries(linewidth=1.5, zorder=1)
        m.drawrivers(color='lightblue', zorder=1)
        return m

    def statistics(self, fn):
        a = np.load(fn)
        lats = a['lat'][:, :, 0]
        lons = a['lon'][:, :, 0]
        times = np.median(a['ttP'], axis=-1)
        tree = spatial.KDTree(zip(lats.ravel(), lons.ravel()))
        vals = []
        perc_max = 84
        perc_min = 16
        for _evid in self.del_coord.keys():
            lon, lat, dt, dm = self.del_coord[_evid]
            distance, index = tree.query(np.array([[lat, lon]]))
            irow, icol = divmod(index[0], lats.shape[1])
            med = np.median(times[:, irow, icol])
            lb = scoreatpercentile(times[:, irow, icol], perc_min)
            ub = scoreatpercentile(times[:, irow, icol], perc_max)
            val = (dt - lb) / (ub - lb)
            if False:
                plt.figure()
                plt.hist(times[:, irow, icol], bins=np.arange(0, 30), normed=True, histtype='step')
                plt.show()
            vals.append(val)

        for _evid in self.del_coord_rt.keys():
            lon, lat, dt, dm = self.del_coord_rt[_evid]
            distance, index = tree.query(np.array([[lat, lon]]))
            irow, icol = divmod(index[0], lats.shape[1])
            med = np.median(times[:, irow, icol])
            lb = scoreatpercentile(times[:, irow, icol], perc_min)
            ub = scoreatpercentile(times[:, irow, icol], perc_max)
            val = (dt - lb) / (ub - lb)
            vals.append(val)
            if False:
                plt.figure()
                plt.hist(times[:, irow, icol, :], bins=np.arange(0, 30), normed=True, histtype='step')
                plt.show()
        idx = np.where((np.array(vals) <= 1.0) & (np.array(vals) >= 0.0))
        print "%.1f lie within the %d and %d percentile" % ((float(idx[0].size) / len(vals)) * 100, perc_min, perc_max)
        # plt.hist(vals, bins=20)
        # plt.show()

    def alert_times_map(self, fout=None, m=None, fig=None, ax=None, scale=10000.,
                         cb=False, magerr=False, interactive=False, msscale=1,
                         realtime=False, cmapname='jet'):
        cmap = cm.ScalarMappable(norm=Normalize(vmin=6, vmax=25), cmap=cmapname)
        if m is None and fig is None and ax is None:
            fig = plt.figure()
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            m = self.background_map(ax)
        dataX = []
        dataY = []
        values = []
        for evid in self.del_coord.keys():
            lon, lat, delay, magdiff = self.del_coord[evid]
            x, y = m(lon, lat)
            dataX.append(x)
            dataY.append(y)
            values.append('%s: %.2f' % (evid, delay))
            cl = cmap.to_rgba(delay)
            if magerr:
                factor = math.sqrt(abs(float(magdiff)))
                sl2 = scale * factor
                p2 = Wedge((x, y), sl2, 0, 360, facecolor=cl,
                           edgecolor='black', picker=5, lw=1.0)
                ax.add_patch(p2)
            else:
                m.plot(x, y, ms=8 * msscale, c=cl, marker='s', picker=5.)
        print "Number of playback events plotted: %d" % len(self.del_coord.keys())

        if realtime:
            for evid in self.del_coord_rt.keys():
                lon, lat, delay, magdiff = self.del_coord_rt[evid]
                x, y = m(lon, lat)
                dataX.append(x)
                dataY.append(y)
                values.append('%s: %.2f' % (evid, delay))
                cl = cmap.to_rgba(delay)
                if magerr:
                    factor = math.sqrt(abs(float(magdiff)))
                    sl2 = scale * factor
                    p2 = Wedge((x, y), sl2, 0, 360, facecolor=cl,
                               edgecolor='black', picker=5, lw=1.0)
                    ax.add_patch(p2)
                else:
                    m.plot(x, y, ms=8 * msscale, c=cl, marker='o', picker=5.)
            print "Number of realtime events plotted: %d" % \
            len(self.del_coord_rt.keys())

        if interactive:
            try:
                # add a pop-up window showing the station and its value
                tooltip = wx.ToolTip(tip='')
                tooltip.Enable(False)
                tooltip.SetDelay(0)
                fig.canvas.SetToolTip(tooltip)
                def onMotion(event):
                    dir(event.artist)
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

        if magerr:
            # Plot legend
            yshift = 23000
            fontsize = 12
            x, y = m(7.8, 45.1)
            ax.text(x, y + yshift, '0.2', horizontalalignment='center',
                    verticalalignment='center', fontsize=fontsize, fontweight='bold')
            p = Wedge((x, y), scale * math.sqrt(0.2), 0, 360, facecolor='black',
                       edgecolor='black', lw=1.0)
            ax.add_patch(p)
            x, y = m(8.6, 45.1)
            ax.text(x, y + yshift, '0.5', horizontalalignment='center',
                    verticalalignment='center', fontsize=fontsize, fontweight='bold')
            p = Wedge((x, y), scale * math.sqrt(0.5), 0, 360, facecolor='black',
                       edgecolor='black', lw=1.0)
            ax.add_patch(p)
            x, y = m(9.3, 45.1)
            ax.text(x, y + yshift, '1.0', horizontalalignment='center',
                    verticalalignment='center', fontsize=fontsize, fontweight='bold')
            p = Wedge((x, y), scale * 1.0, 0, 360, facecolor='black',
                       edgecolor='black', lw=1.0)
            ax.add_patch(p)
            ax.text(x + 200000, y , 'Magnitude error', horizontalalignment='right',
                    verticalalignment='center', fontsize=fontsize, fontweight='bold')

        if cb:
            # Colorbar
            cax = fig.add_axes([0.87, 0.1, 0.05, 0.8])
            cb = ColorbarBase(cax, cmap=cmapname,
                              norm=Normalize(vmin=6., vmax=25.))
            cb.set_label('Time since origin time [s]')
        if fout is not None:
            fig.savefig(fout, dpi=300, bbox_inches='tight')
