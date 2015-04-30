#!/usr/bin/env python
"""
Compute theoretical alert delays.

Created on Jan 22, 2015

@author: behry
"""
import os
import sys
sys.path.append(os.path.join(os.environ['HOME'], 'mypy'))
import progressbar as pg
import numpy as np
import pyproj
import ipdb


class NetworkInfo:
    """
    Dummy implementation for a class that returns station names
    and coordinates, and network codes.
    """
    def __init__(self):
        self.lat = []
        self.lon = []
        self.nm = []
        self.nw = []
        self.size = 0

    def read(self):
        pass

    def get_lat(self):
        pass

    def get_lon(self):
        pass

    def get_names(self):
        pass

    def get_networkcodes(self):
        pass


class DelayObs:
    """
    Dummy class to load processing latencies.
    """
    def __init__(self):
        pass

    def get_pick_delays(self):
        pass

    def get_pick_default(self):
        pass

    def get_envelope_default(self):
        pass

    def get_envelope_delays(self):
        pass

    def get_associator_delays(self):
        pass

    def get_magnitude_delays(self):
        pass


class DelayException(Exception):
    pass


class DelayEEW:

    def __init__(self):
        pass

    def compute(self, networkinfo, elon, elat, edep, vp=6.5, vs=3.5, nnst=6,
                procdelay=False, target=None, nmaps=500,
                resultsfn='p_wave_tt_6_stations.npz', datadir='./data',
                latencies=None, vp3d=False):

        if target is not None:
            tln, tlt = target

        networks = networkinfo.get_networks()

        if procdelay and latencies is not None:
            pkdel = latencies.get_pick_delays()
            envdel = latencies.get_envelope_delays()
            ascdel = latencies.get_associator_delays()
            magdel = latencies.get_magnitude_delays()
            pkdefault = latencies.get_pick_default()
            envdefault = latencies.get_envelope_default()

        # do we have a grid or a set of points?
        if elon.shape != elat.shape or elon.shape != edep.shape:
            raise DelayException("elon, elat, and edep must have the same shape!")

        outshape = [nmaps]
        for _s in elon.shape:
            outshape.append(_s)

        lat = elat.ravel()
        lon = elon.ravel()
        dep = edep.ravel()
        ngp = lon.size

        ttP = np.zeros(outshape)
        ttPtmp = np.zeros(ngp)
        tstarget = ttPtmp.copy()

        g = pyproj.Geod(ellps='sphere')
        # Setup progressbar
        widgets = ['tt: ', pg.Percentage(), ' ', pg.Bar('#'),
                   ' ', pg.ETA()]
        pbar = pg.ProgressBar(widgets=widgets, maxval=nmaps * ngp).start()
        no_env_dl = []
        no_pk_dl = []
        for nm in range(nmaps):
            idx = 0
            for _lat, _lon, _dep in zip(lat, lon, dep):
                pbar.update(nm * ngp + idx)
                min_dt = 1.e38
                for net in networks.keys():
                    if len(networks[net]['lat']) < 1:
                        continue
                    stlat = networks[net]['lat']
                    stlon = networks[net]['lon']
                    nwcode = np.array(networks[net]['nw'])
                    names = np.array(networks[net]['nm'])
                    maxnstat = min(10, len(stlat))
                    nnst = min(nnst, len(stlat))

                    if vp3d:
                        datadir = './data/ttime_ch/'
                        fin = os.path.join(datadir, "%.4f_%.4f_%.1f.npy" % \
                                            (_lat, _lon, _dep))
                        if not os.path.isfile(fin):
                            raise DelayException('File %s does not exist!' % fin)
                        a = np.load(fin)
                        pt = a['ttime'][:, 0]
                        dt = max(pt[0:nnst])
                        stat_names = a['name']
                    else:
                        # Find the <nnst> nearest stations
                        lats = np.ones((len(stlat),)) * _lat
                        lons = np.ones((len(stlon),)) * _lon
                        az, baz, dist = g.inv(lons, lats, stlon, stlat)
                        dist_sorted = np.sort(dist)
                        stat_names = names[np.argsort(dist)]
                        dz = np.ones((maxnstat,)) * _dep
                        pt = np.sqrt(dz * dz + dist_sorted[0:maxnstat] / 1000. * dist_sorted[0:maxnstat] / 1000.) / vp
                        dt = max(pt[0:nnst])
                    if target is not None:
                        azt, bazt, distt = g.inv(_lon, _lat, tln, tlt)
                        distt /= 1000.
                        tstarget[idx] = np.sqrt(distt * distt + _dep * _dep) / vs

                    # Add delays
                    if procdelay:
                        pk_delays = []
                        env_delays = []
                        for stat in stat_names[0:maxnstat]:
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

                        if len(pk_delays) != maxnstat or len(env_delays) != maxnstat:
                            raise DelayException('Number of delays is not equal to %d' % nnst)
                        temp = np.sort(pt + np.array(pk_delays)) + ascdel[np.random.randint(0, len(ascdel))]
                        origin_delay = max(temp[0:nnst])
                        waveform_delay = min(pt + np.array(env_delays))
                        dt = max(origin_delay, waveform_delay) + magdel[np.random.randint(0, len(magdel))]
                    if dt < min_dt:
                        min_dt = dt
                ttPtmp[idx] = min_dt
                idx += 1
            ttP[nm] = ttPtmp.reshape(elon.shape)

        tstarget = tstarget.reshape(elon.shape)
        if resultsfn is not None:
            np.savez(resultsfn, ttP=ttP, tstarget=tstarget, lat=elat, lon=elon,
                     dep=edep)
        pbar.finish()
        print "No envelope delay data available for the following stations:"
        print ' '.join(no_env_dl)
        print "No pick delay info available for the following stations:"
        print ' '.join(no_pk_dl)
        return elat, elon, edep, ttP, tstarget

if __name__ == '__main__':
    pass
