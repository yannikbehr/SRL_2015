#!/usr/bin/env python
"""
Created on Oct 18, 2013

@author: behry
"""

import math
import sys
sys.path.append('/home/behry/workspace/eew/reports')
from obspy import UTCDateTime
import matplotlib
try:
    matplotlib.use('WXAgg')
    import wx
except:
    print "WX package for Python not installed"
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import scoreatpercentile
from mpl_toolkits.basemap import Basemap
from scipy.io.netcdf import netcdf_file as Dataset #from Scientific.IO.NetCDF import NetCDFFile as Dataset
from matplotlib.colors import LightSource
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.pyplot import cm
from matplotlib.patches import Wedge
from reports_parser import ReportsParser
from obspy.core.util import gps2DistAzimuth
from point_in_polygon import EventCA, EventSoCal
from scipy import spatial
import pyproj
from delayeew import DelayEEW


class AlertTimes:
    """
    Analyse observed and predicted alert times in California.
    """

    def __init__(self):
        self.del_coord = []
        self.event_excludes = ['NC_71736656']

    def background_map(self, ax):
        llcrnrlat, urcrnrlat, llcrnrlon, urcrnrlon, lat_ts = (31, 44, -126, -113, 37.5)
        m = Basemap(projection='merc', llcrnrlat=llcrnrlat,
                    urcrnrlat=urcrnrlat, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                    lat_ts=lat_ts, resolution='i', ax=ax)
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
        m.drawstates()
        m.drawrivers(color='lightblue', zorder=1)
        return m

    def popup(self, fig, dataX, dataY, values):
        # add a pop-up window showing the station and its value
        try:
            tooltip = wx.ToolTip(tip='')
            tooltip.Enable(False)
            tooltip.SetDelay(0)
            fig.canvas.SetToolTip(tooltip)

            def onMotion(event):
                # dir(event.artist)
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

    def load_event_info(self, fn):
        ei = {}
        fh = open(fn)
        for _l in fh.readlines():
            a = _l.split()
            # if int(a[1]) != 0
            ei[a[0]] = int(a[1])
        return ei

    def closest_stations(self, stations, lat, lon, nnst=4, vp=6.5, depth=8):
        """
        Find the 4 closest stations that are most likely to have detected this
        event first.
        """
        g = pyproj.Geod(ellps='WGS84')
        min_dt = 1.e38
        for net in ['bk', 'ca', 'mp']:
            if len(stations.networks[net]['lat']) < 1:
                continue
            stlat = stations.networks[net]['lat']
            stlon = stations.networks[net]['lon']
            nwcode = np.array(stations.networks[net]['nw'])
            names = np.array(stations.networks[net]['nm'])
            # Find the <nnst> nearest stations
            lats = np.ones((len(stlat),)) * lat
            lons = np.ones((len(stlon),)) * lon
            az, baz, dist = g.inv(lons, lats, stlon, stlat)
            dist_sorted = np.sort(dist)
            stat_names = names[np.argsort(dist)]
            networks = nwcode[np.argsort(dist)]
            dz = np.ones((10,)) * depth
            distance = np.sqrt(dz * dz + dist_sorted[0:10] / 1000. * dist_sorted[0:10] / 1000.)
            dt = max(distance[0:nnst] / vp)
            if dt < min_dt:
                min_dt = dt
                # stats = ['%s.%s' % (net, st) for st in stat_names[0:nnst]]
                stats = stat_names[0:nnst]
        return stats

    def statistics(self, fns, fn, stationfn, eventinfo=None, latencies=None,
                   computedelay=False, map=False, interactive=False):
        """
        Compare predicted and observed alert times quantitatively.
        """
        a = np.load(fn)
        lats_tt = a['lat'][:, :, 0]
        lons_tt = a['lon'][:, :, 0]
        times = np.median(a['ttP'], axis=-1)
        tree = spatial.KDTree(zip(lats_tt.ravel(), lons_tt.ravel()))
        vals = []
        perc_max = 84
        perc_min = 16

        rp = ReportsParser(dmin=UTCDateTime(2012, 1, 1, 0, 0, 0),
                           dmax=UTCDateTime(2013, 11, 1, 0, 0, 0))
        # t = EventCA()
        t = EventSoCal()
        rp.sfilter = t.point_in_polygon

        for _f in fns:
            rp.read_reports(_f)

        correct = rp.get_correct(mmin=3.5, mmax=10.0)
        pid = correct[:, 0]
        ot = correct[:, 2].astype('float')
        lats = correct[:, 3].astype('float')
        lons = correct[:, 4].astype('float')
        deps = correct[:, 5].astype('float')
        mags = correct[:, 6].astype('float')
        ts1 = correct[:, 7].astype('float')
        lats1 = correct[:, 9].astype('float')
        lons1 = correct[:, 10].astype('float')
        mags1 = correct[:, 12].astype('float')
        rfns = correct[:, 21]
        diff = ts1 - ot
        magdiff = mags - mags1
        cnt = 0
        allcnt = 0
        allm = []
        dataX = []
        dataY = []
        popup_values = []

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        m = self.background_map(ax)
        cmap = cm.ScalarMappable(norm=Normalize(vmin=0, vmax=2), cmap='RdBu_r')
        stats_used = []
        for lon, lat, dep, delay, evid, lat1, lon1, dmag, time, mag, rfn in \
            zip(lons, lats, deps, diff, pid, lats1, lons1, magdiff, ot, mags, rfns):
            allcnt += 1
            try:
                if eventinfo is not None and len(eventinfo[evid]) != 4:
                    # print "Event %s does not have 4 initial picks." % evid
                    continue
            except KeyError:
                print "No event information available for: %s (%s)" % (evid, UTCDateTime(time))
                continue
            if evid in self.event_excludes:
                print "Event %s was set to be excluded." % evid
                continue

            if computedelay:
                # Compute the expected alert time for the actual epicenter and
                # the first stations that detected the event
                class NetworkInfo:
                    def __init__(self):
                        self.networks = {'ca':{'lat': [], 'lon': [], 'chn': [],
                                               'nw': [], 'nm': [], 'lc': [],
                                               'color':'black',
                                               'label':'UC Berkeley'}}
                    def get_networks(self):
                        return self.networks

                # read in SCEDC master station list
                fh = open(stationfn)
                scedc_stations = {}
                for _l in fh.readlines():
                    if _l.startswith('#'):
                        continue
                    net, sta, chan, loc, lt, ln, elev, ondate, offdate = _l.split()
                    ns = '.'.join((net, sta))
                    if ns not in scedc_stations:
                        scedc_stations[ns] = (float(lt), float(ln))
                ni = NetworkInfo()
                for _st in eventinfo[evid]:
                    ni.networks['ca']['lat'].append(scedc_stations[_st][0])
                    ni.networks['ca']['lon'].append(scedc_stations[_st][1])
                    ni.networks['ca']['nm'].append(_st)
                    if _st not in stats_used:
                        stats_used.append(_st)
                de = DelayEEW()
                elat, elon, edep, ttP, tstarget = \
                de.compute(ni, np.array([float(lon)]), np.array([float(lat)]),
                           np.array([float(dep)]),
                           vp=6.5, vs=3.5, nnst=4, procdelay=True, nmaps=500,
                           resultsfn=None, latencies=latencies)
                med = np.median(ttP)
                lb = scoreatpercentile(ttP, perc_min)
                ub = scoreatpercentile(ttP, perc_max)
            else:
                distance, index = tree.query(np.array([[lat, lon]]))
                irow, icol = divmod(index[0], lats_tt.shape[1])
                med = np.median(times[:, irow, icol])
                lb = scoreatpercentile(times[:, irow, icol], perc_min)
                ub = scoreatpercentile(times[:, irow, icol], perc_max)

            cnt += 1
            allm.append(mag)
            val = (delay - lb) / (ub - lb)
            print med, lb, ub, delay, val, med - delay
            vals.append(val)
            cl = cmap.to_rgba(val)
            x, y = m(lon, lat)
            dataX.append(x)
            dataY.append(y)
            info = '%s: %.2f %s\n' % (UTCDateTime(time), mag, evid)
            info += '%.2f %.2f %.2f\n' % (delay, med, val)
            for _st in eventinfo[evid]:
                info += ' %s' % _st
            popup_values.append(info)
            m.plot(x, y, ms=8, c=cl, marker='o', picker=5.)
            # plt.figure()
            # plt.hist(times[ilon, ilat, :], bins=np.arange(0, 30), normed=True, histtype='step')
            # plt.show()
        print "Stations used in detections:"
        print stats_used
        idx = np.where((np.array(vals) <= 1.0) & (np.array(vals) >= 0))
        print "%.1f lie within the %d and %d percentile" % ((idx[0].size / float(len(vals))) * 100, perc_min, perc_max)
        # plt.plot(allm, vals, 'bo')
        if interactive:
            self.popup(fig, dataX, dataY, popup_values)
        cax = fig.add_axes([0.87, 0.1, 0.05, 0.8])
        cb = ColorbarBase(cax, cmap='RdBu_r',
                          norm=Normalize(vmin=0, vmax=2))
        cb.set_label('Alert accuracy')

        plt.figure()
        plt.hist(vals, bins=20)
        plt.show()

    def alert_times_map(self, fns, m=None, fig=None, ax=None, scale=10000.,
                        cb=True, disterr=False, interactive=False,
                        eventinfo=None, msscale=1, cmapname='jet'):
        """
        Plot a map of observed alert times.
        """
        cmap = cm.ScalarMappable(norm=Normalize(vmin=6, vmax=25), cmap=cmapname)
        rp = ReportsParser(dmin=UTCDateTime(2012, 1, 1, 0, 0, 0),
                           dmax=UTCDateTime(2013, 11, 1, 0, 0, 0))
        t = EventCA()
        rp.sfilter = t.point_in_polygon

        for _f in fns:
            rp.read_reports(_f)

        correct = rp.get_correct(mmin=3.5, mmax=10.0)
        pid = correct[:, 0]
        ot = correct[:, 2].astype('float')
        lats = correct[:, 3].astype('float')
        lons = correct[:, 4].astype('float')
        mags = correct[:, 6].astype('float')
        ts1 = correct[:, 7].astype('float')
        lats1 = correct[:, 9].astype('float')
        lons1 = correct[:, 10].astype('float')
        mags1 = correct[:, 12].astype('float')
        rfns = correct[:, 21]
        diff = ts1 - ot
        magdiff = mags - mags1

        if m is None and fig is None and ax is None:
            fig = plt.figure()
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
            m = self.background_map(ax)
        dataX = []
        dataY = []
        values = []
        # load event info
        cnt = 0
        allcnt = 0
        for lon, lat, delay, evid, lat1, lon1, dmag, time, mag, rfn in \
            zip(lons, lats, diff, pid, lats1, lons1, magdiff, ot, mags, rfns):
            allcnt += 1
            try:
                if eventinfo is not None and len(eventinfo[evid]) != 4:
                    print "Event %s does not have 4 initial picks." % evid
                    continue
            except KeyError:
                print "No event information available for: %s (%s)" % (evid, UTCDateTime(time))
                continue
            if evid in self.event_excludes:
                print "Event %s was set to be excluded." % evid
                continue
            cnt += 1
            ddist, az, baz = gps2DistAzimuth(lat, lon, lat1, lon1)
            ddist /= 1000.
            x, y = m(lon, lat)
            dataX.append(x)
            dataY.append(y)
            info = '%s: %.2f %.2f %s' % (UTCDateTime(time), delay, mag, evid)
            for _st in eventinfo[evid]:
                info += ' %s' % _st
            values.append(info)
            cl = cmap.to_rgba(delay)
            if disterr:
                factor = math.sqrt(abs(float(ddist)))
                sl2 = scale * factor
                p2 = Wedge((x, y), sl2, 0, 360, facecolor=cl,
                           edgecolor='black', picker=5, lw=1.0)
                ax.add_patch(p2)
            else:
                m.plot(x, y, ms=8 * msscale, c=cl, marker='o', picker=5.)
        print "Plotted %d out of %d events." % (cnt, allcnt)
        if interactive:
            self.popup(fig, dataX, dataY, values)
        if cb:
            # Colorbar
            cax = fig.add_axes([0.87, 0.1, 0.05, 0.8])
            cb = ColorbarBase(cax, cmap=cmapname,
                              norm=Normalize(vmin=6., vmax=25.))
            cb.set_label('Time since origin time [s]')


if __name__ == '__main__':
    pass
