#!/usr/bin/env python
"""
Make a plot of P-wave travel time to 6 stations in the network.

Created on Jun 24, 2013

@author: behry
"""
import os
import json
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'mypy'))
import progressbar as pg
import numpy as np
import pyproj


class StationInfo:

    def __init__(self):
        self.networks = {}
        self.networks['bk'] = {'lat': [], 'lon': [], 'chn': [],
                               'nw': [], 'nm': [], 'lc': []}
        self.networks['ca'] = {'lat': [], 'lon': [], 'chn': [],
                               'nw': [], 'nm': [], 'lc': []}
        self.networks['mp'] = {'lat': [], 'lon': [], 'chn': [],
                               'nw': [], 'nm': [], 'lc': []}
        self.lat = []
        self.lon = []
        self.chn = []
        self.nw = []
        self.nm = []
        self.lc = []
        self.excludes = ['TSCN', 'SRI', 'OXMT', 'MHDL', 'KRP', 'SCI2', 'OSI',
                         'GR2', 'MSC', 'DEV', 'MONP2', 'MGE', 'MTP', 'ERR',
                         'DRE', 'MNRC', 'KRMB', 'KSXB', 'CTD', 'KCC', 'MMLB',
                         'PASC', 'LUG']

    def read(self, fn, sm=True, statewide=False):
        f = open(fn, 'r')
        for line in f.readlines():
            nt, st, chn, loc, lat, lon = line.split()
            if not sm:
                if chn[0:2] == 'HG' or chn[0:2] == 'HN':
                    continue
            if st in self.excludes:
                continue
            ns = '%s.%s' % (nt, st)
            # To simulate a statewide system all stations are assigned
            # to a single network
            if statewide:
                if ns not in self.networks['ca']['nm']:
                    self.networks['ca']['nm'].append(ns)
                    self.networks['ca']['nw'].append(nt)
                    self.networks['ca']['chn'].append(chn)
                    self.networks['ca']['lc'].append(loc)
                    self.networks['ca']['lat'].append(float(lat))
                    self.networks['ca']['lon'].append(float(lon))
            else:
                if nt == 'AZ' or nt == 'CI':
                    if ns not in self.networks['ca']['nm']:
                        self.networks['ca']['nm'].append(ns)
                        self.networks['ca']['nw'].append(nt)
                        self.networks['ca']['chn'].append(chn)
                        self.networks['ca']['lc'].append(loc)
                        self.networks['ca']['lat'].append(float(lat))
                        self.networks['ca']['lon'].append(float(lon))
                elif nt == 'NP' or nt == 'NC':
                    if ns not in self.networks['mp']['nm']:
                        self.networks['mp']['nm'].append(ns)
                        self.networks['mp']['nw'].append(nt)
                        self.networks['mp']['chn'].append(chn)
                        self.networks['mp']['lc'].append(loc)
                        self.networks['mp']['lat'].append(float(lat))
                        self.networks['mp']['lon'].append(float(lon))
                elif nt == 'BK':
                    if ns not in self.networks['bk']['nm']:
                        self.networks['bk']['nm'].append(ns)
                        self.networks['bk']['nw'].append(nt)
                        self.networks['bk']['chn'].append(chn)
                        self.networks['bk']['lc'].append(loc)
                        self.networks['bk']['lat'].append(float(lat))
                        self.networks['bk']['lon'].append(float(lon))
                else:
                    print "Unknown network: %s %s " % (nt, st)
        f.close()


def load_latencies(datadir, perstation=False):
    # get all the delay distributions
    pkdelfn = os.path.join(datadir, 'single_station_pk_delays_ca.txt')
    envdelfn = os.path.join(datadir, 'single_station_env_delays_ca.txt')
    octdelfn = os.path.join(datadir, 'eewvs_origin_delays.npz')
    magdelfn = os.path.join(datadir, 'magnitude_estimation_delay.npz')

    fh = open(pkdelfn)
    pkdel = json.load(fh)
    fh.close()
    bkdelays = []
    mpdelays = []
    cadelays = []
    for _k in pkdel.keys():
        net, stat = _k.split('.')
        if net == 'NP' or net == 'NC':
            if perstation:
                mpdelays.append(np.median(pkdel[_k]))
            else:
                mpdelays += pkdel[_k]
        elif net == 'CI' or net == 'AZ':
            if perstation:
                cadelays.append(np.median(pkdel[_k]))
            else:
                cadelays += pkdel[_k]
        elif net == 'BK':
            if perstation:
                bkdelays.append(np.median(pkdel[_k]))
            else:
                bkdelays += pkdel[_k]
        else:
            print "Unknown network %s" % net

    pkdefault = {'mp': np.median(mpdelays),
                 'bk': np.median(bkdelays),
                 'ca': np.median(cadelays)}

    fh = open(envdelfn)
    envdel = json.load(fh)
    fh.close()
    bkdelays = []
    mpdelays = []
    cadelays = []
    for _k in envdel.keys():
        net, stat = _k.split('.')
        if net == 'NP' or net == 'NC':
            if perstation:
                mpdelays.append(np.median(envdel[_k]))
            else:
                mpdelays += envdel[_k]
        elif net == 'CI' or net == 'AZ':
            if perstation:
                cadelays.append(np.median(envdel[_k]))
            else:
                cadelays += envdel[_k]
        elif net == 'BK':
            if perstation:
                bkdelays.append(np.median(envdel[_k]))
            else:
                bkdelays += envdel[_k]
        else:
            print "Unknown network %s" % net

    envdefault = {'mp': np.median(mpdelays),
                 'bk': np.median(bkdelays),
                 'ca': np.median(cadelays)}

    a = np.load(octdelfn)
    octdel = a['acdl']
    a = np.load(magdelfn)
    magdel = a['delays']

    return pkdel, envdel, octdel, magdel, pkdefault, envdefault


def traveltime(stationlist, depth=8, vp=6.5, vs=3.5, nnst=6, procdelay=False,
               sm=False, boxin=(32, 43, -125, -114), delayonly=False,
               nmaps=500, resultsfn='p_wave_tt_4_stations_ca.npz',
               datadir='./data', statewide=False):

    latmin, latmax, lonmin, lonmax = boxin
    # Number of gridpoints in either direction
    ngp = 201
    lat = np.linspace(latmin, latmax, ngp)
    lon = np.linspace(lonmin, lonmax, ngp)

    # Station list
    stFileName = stationlist
    station = StationInfo()
    station.read(stFileName, sm=sm, statewide=statewide)

    if procdelay:
        pkdel, envdel, octdel, magdel, pkdefault, envdefault = \
        load_latencies(datadir)

    # Loop over all gridpoints
    ttP = np.zeros((ngp, ngp, nmaps))

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
            for _nm in xrange(nmaps):
                min_dt = 1.e38
                for net in ['bk', 'ca', 'mp']:
                    if len(station.networks[net]['lat']) < 1:
                        continue
                    stlat = station.networks[net]['lat']
                    stlon = station.networks[net]['lon']
                    nwcode = np.array(station.networks[net]['nw'])
                    names = np.array(station.networks[net]['nm'])
                    # Find the <nnst> nearest stations
                    lats = np.ones((len(stlat),)) * lat[ilat]
                    lons = np.ones((len(stlon),)) * lon[ilon]
                    az, baz, dist = g.inv(lons, lats, stlon, stlat)
                    dist_sorted = np.sort(dist)
                    stat_names = names[np.argsort(dist)]
                    networks = nwcode[np.argsort(dist)]
                    dz = np.ones((10,)) * depth
                    distance = np.sqrt(dz * dz + dist_sorted[0:10] / 1000. * dist_sorted[0:10] / 1000.)
                    dt = max(distance[0:nnst] / vp)

                    # Add delays
                    if procdelay:
                        pk_delays = []
                        env_delays = []
                        for stat in stat_names[0:10]:
                            if stat in pkdel.keys():
                                pk_delays.append(pkdel[stat][np.random.randint(0, len(pkdel[stat]))])
                            else:
                                if stat not in no_pk_dl:
                                    no_pk_dl.append(stat)
                                pk_delays.append(pkdefault[net])

                            if stat in envdel.keys():
                                env_delays.append(envdel[stat][np.random.randint(0, len(envdel[stat]))])
                            else:
                                if stat not in no_env_dl:
                                    no_env_dl.append(stat)
                                env_delays.append(envdefault[net])

                        if len(pk_delays) != 10 or len(env_delays) != 10:
                            raise Exception('Number of delays is not equal to %d (%d,%d)' % (nnst, len(pk_delays), len(env_delays)))
                        temp = np.sort(distance / vp + np.array(pk_delays)) + octdel[np.random.randint(0, len(octdel))]
                        origin_delay = max(temp[0:nnst])
                        waveform_delay = min(distance / vp + np.array(env_delays))
                        if delayonly:
                            origin_delay = max(np.array(pk_delays[0:nnst])) + magdel[np.random.randint(0, len(magdel))]
                            waveform_delay = min(np.array(env_delays))
                        dt = max(origin_delay, waveform_delay) + magdel[np.random.randint(0, len(magdel))]
                        # print "origin delay: %.2f; waveform delay: %.2f; delay time: %.2f" % (origin_delay, waveform_delay, dt)

                    if dt < min_dt:
                        min_dt = dt

                ttP[ilon, ilat, _nm] = min_dt

    np.savez(resultsfn, ttP=ttP, lat=lat, lon=lon)
    pbar.finish()
    print "No envelope delay data available for the following stations:"
    print ' '.join(no_env_dl)
    print "No picker delay info available for the following stations:"
    print ' '.join(no_pk_dl)

if __name__ == '__main__':
    stationlist = 'data/stations_ca.txt'
    resultsfn = 'data/p_wave_tt_4_stations_ca_procdel.npz'
    traveltime(stationlist, procdelay=True, nmaps=500, new=True,
                  resultsfn=resultsfn)
