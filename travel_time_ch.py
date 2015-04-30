#!/usr/bin/env python
"""
Make a plot of P-wave travel time to 6 stations in the network.

Created on Jun 24, 2013

@author: behry
"""
import os
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'mypy'))
import progressbar as pg
import numpy as np
import pyproj
import json


class StationInfo:

    def __init__(self):
        self.lat = []
        self.lon = []
        self.chn = []
        self.nw = []
        self.nm = []
        self.lc = []
        self.excludes = ['RA.OGSI', 'RA.STBO', 'RA.STFL']

    def read(self, fn, sm=True):
        f = open(fn, 'r')
        for line in f.readlines():
            nt, st, chn, loc, lat, lon = line.split()
            if not sm:
                if chn[0:2] == 'HG' or chn[0:2] == 'HN':
                    continue
            ns = '%s.%s' % (nt, st)
            if ns in self.excludes:
                continue
            if ns not in self.nm:
                self.nm.append(ns)
                self.nw.append(nt)
                self.chn.append(chn)
                self.lc.append(loc)
                self.lat.append(float(lat))
                self.lon.append(float(lon))
        f.close()


def load_latencies(datadir, perstation=False):
    # get all the delay distributions
    pkdelfn = os.path.join(datadir, 'single_station_pk_delays_ch.txt')
    envdelfn = os.path.join(datadir, 'single_station_env_delays_ch.txt')
    octdelfn = os.path.join(datadir, 'origin_estimation_delay_ch.npz')
    magdelfn = os.path.join(datadir, 'magnitude_estimation_delay.npz')

    fh = open(pkdelfn)
    pkdel = json.load(fh)
    fh.close()
    delays = []
    for st in pkdel.keys():
        if perstation:
            delays.append(np.median(pkdel[st]))
        else:
            delays += pkdel[st]
    pkdefault = np.median(delays)

    fh = open(envdelfn)
    envdel = json.load(fh)
    fh.close()
    delays = []
    for st in envdel.keys():
        if perstation:
            delays.append(np.median(envdel[st]))
        else:
            delays += envdel[st]
    envdefault = np.median(delays)

    a = np.load(octdelfn)
    octdel = a['delays']
    a = np.load(magdelfn)
    magdel = a['delays']

    return pkdel, envdel, octdel, magdel, pkdefault, envdefault


def traveltime(stationlist, depth=8, vp=6.5, vs=3.5, nnst=6, procdelay=False,
               sm=False, boxin=(45.4, 48.3, 5.6, 11.1), target=None, nmaps=500,
               new=False, resultsfn='p_wave_tt_6_stations.npz',
               datadir='./data'):

    if target is not None:
        tln, tlt = target

    latmin, latmax, lonmin, lonmax = boxin
    # Number of gridpoints in either direction
    ngp = 201
    lat = np.linspace(latmin, latmax, ngp)
    lon = np.linspace(lonmin, lonmax, ngp)

    # Station list
    stFileName = stationlist
    station = StationInfo()
    station.read(stFileName, sm=sm)

    stlat = station.lat
    stlon = station.lon
    names = np.array(station.nm)

    if procdelay:
        pkdel, envdel, octdel, magdel, pkdefault, envdefault = \
        load_latencies(datadir)

    # Loop over all gridpoints
    ttP = np.zeros((ngp, ngp, nmaps))
    tstarget = np.zeros((ngp, ngp, nmaps))

    g = pyproj.Geod(ellps='sphere')
    # Setup progressbar
    widgets = ['tt: ', pg.Percentage(), ' ', pg.Bar('#'),
               ' ', pg.ETA()]
    pbar = pg.ProgressBar(widgets=widgets, maxval=ngp).start()

    no_env_dl = []
    no_pk_dl = []
    for ilon in range(ngp):
        pbar.update(ilon)
        for ilat in range(ngp):
            # Find the <nnst> nearest stations
            lats = np.ones((len(stlat),)) * lat[ilat]
            lons = np.ones((len(stlon),)) * lon[ilon]
            az, baz, dist = g.inv(lons, lats, stlon, stlat)
            dist_sorted = np.sort(dist)
            stat_names = names[np.argsort(dist)]
            dz = np.ones((10,)) * depth
            distance = np.sqrt(dz * dz + dist_sorted[0:10] / 1000. * dist_sorted[0:10] / 1000.)
            dt = max(distance[0:nnst] / vp)
            if target is not None:
                azt, bazt, distt = g.inv(lon[ilon], lat[ilat], tln, tlt)
                distt /= 1000.
                tstarget[ilon, ilat, 0] = np.sqrt(distt * distt + depth * depth) / vs

            # Add delays
            if procdelay:
                for _nm in xrange(nmaps):
                    pk_delays = []
                    env_delays = []
                    for stat in stat_names[0:10]:
                        if True:
                            pk_delays.append(np.random.uniform(0.5, 2.5))
                            env_delays.append(np.random.uniform(0.5, 2.5))
                        if False:
                            if stat in pkdel.keys():
                                pk_delays.append(pkdel[stat][np.random.randint(0, len(pkdel[stat]))])
                            else:
                                if stat not in no_pk_dl:
                                    no_pk_dl.append(stat)
                                pk_delays.append(pkdefault)

                            if stat in envdel.keys():
                                env_delays.append(envdel[stat][np.random.randint(0, len(envdel[stat]))])
                            else:
                                if stat not in no_env_dl:
                                    no_env_dl.append(stat)
                                env_delays.append(envdefault)

                    if len(pk_delays) != 10 or len(env_delays) != 10:
                        raise Exception('Number of delays is not equal to %d' % nnst)
                    temp = np.sort(distance / vp + np.array(pk_delays)) + octdel[np.random.randint(0, len(octdel))]
                    origin_delay = max(temp[0:nnst])
                    waveform_delay = min(distance / vp + np.array(env_delays))
                    dt = max(origin_delay, waveform_delay) + magdel[np.random.randint(0, len(magdel))]
                    ttP[ilon, ilat, _nm] = dt
            else:
                ttP[ilon, ilat, 0] = dt

    np.savez(resultsfn, ttP=ttP, tstarget=tstarget, lat=lat, lon=lon)
    pbar.finish()
    print "No envelope delay data available for the following stations:"
    print ' '.join(no_env_dl)
    print "No pick delay info available for the following stations:"
    print ' '.join(no_pk_dl)


if __name__ == '__main__':
    pass
