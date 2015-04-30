#!/usr/bin/env python
"""
Re-implementation of the Swiss GMPE based on ??
following the matlab implementation of Carlo Cauzzi.

Created on Feb 18, 2014

@author: behry
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import pyproj

class SwissGmpe:

    def __init__(self):
        ord = 1
        coeff_fn_alpine = os.path.join(os.path.dirname(__file__),
                                'data/sd60_a.dat')
        self.c_alp = np.loadtxt(coeff_fn_alpine, usecols=[ord], unpack=True)
        coeff_fn_foreland = os.path.join(os.path.dirname(__file__),
                                'data/sd60_f.dat')
        self.c_fore = np.loadtxt(coeff_fn_foreland, usecols=[ord], unpack=True)
        self.amp = None
        self.g = pyproj.Geod(ellps='sphere')

    def grid_setup(self, lat, lon):
        self.latgrd, self.longrd = np.meshgrid(lat, lon)
        self.lat = lat
        self.lon = lon
        thresh = 0.39 * self.longrd + 44
        self.forelandmask = np.where(self.latgrd > thresh, 1, 0)
        self.alpinemask = np.where(self.latgrd > thresh, 0, 1)

    def get_amplification(self, new=True):
        fn = os.path.join(os.path.dirname(__file__),
                          'data/ampli_pga_grid.csv')
        savefile = os.path.join(os.path.dirname(__file__),
                                'data/amplification_ch.npz')
        if new:
            amp_alpine = []
            amp_foreland = []
            lta, lna, a = np.loadtxt(fn, skiprows=1, usecols=(0, 1, 5),
                                     unpack=True, delimiter=',')
            self.amp = np.ones((self.lat.size, self.lon.size))
            print self.amp.shape
            for _lt, _ln, _a in zip(lta, lna, a):
                idlt = np.argmin(abs(self.lat - _lt))
                idln = np.argmin(abs(self.lon - _ln))
                self.amp[idlt, idln] = _a
                if 0.39 * _ln + 44 < _lt:
                    amp_foreland.append(_a)
                else:
                    amp_alpine.append(_a)
            print np.mean(amp_alpine), np.mean(amp_foreland)
            np.savez(savefile, amp=self.amp)
        else:
            a = np.load(savefile)
            self.amp = a['amp']

    def maxpga(self, mag, eqlat, eqlon, bz):
        # Distinction between alpine and foreland regions
        az, baz, dist = self.g.inv(self.longrd, self.latgrd,
                              np.ones(self.longrd.shape) * eqlon,
                              np.ones(self.longrd.shape) * eqlat)
        dist /= 1000.
        if eqlat > 0.39 * eqlon + 44:
            r, pga = self.get_pga(mag, dist , 'foreland', 'median')
        else:
            r, pga = self.get_pga(mag, dist / 1000., 'alpine', 'median')
        if self.amp is not None:
            pga *= self.amp
        pga = np.where(dist > bz, pga, 0.)
        return pga.max()


    def get_pga(self, magnitude, distance, region, mode, sn=False):
        # foreland means lat > 0.39*lonn + 44
        # M is MW
        if region == 'alpine':
            c = self.c_alp
        elif region == 'foreland':
            c = self.c_fore
        else:
            print "Region has to be either 'alpine' or 'foreland'!"

        a1 = c[0]
        a2 = c[1]
        a3 = c[2]
        a4 = c[3]
        a5 = c[4]
        a6 = c[5]
        a7 = c[6]
        a8 = c[7]
        a9 = c[8]
        a10 = c[9]
        a11 = c[10]
        a12 = c[11]
        a13 = c[12]
        a14 = c[13]
        a15 = c[14]
        a16 = c[15]
        a17 = c[16]
        a18 = c[17]
        a19 = c[18]
        a20 = c[19]
        a21 = c[20]
        a22 = c[21]
        a23 = c[22]

        M = magnitude
        R = distance

        if region == 'alpine':
            if M > 5:
                rmin = 0.55
            elif M > 4.7:
                rmin = -2.8 * M + 14.55
            else:
                rmin = -0.295 * M + 2.65
        elif region == 'foreland':
            if M > 5.5:
                rmin = 0.55
            elif M > 4.7:
                rmin = -2.067 * M + 11.92
            else:
                rmin = -0.291 * M + 3.48
        else:
            print "Region has to be either 'alpine' or 'foreland'!"

        d = np.log10(np.where(R > rmin, R, rmin))
        y = np.zeros(R.shape)
        ycor = np.zeros(R.shape)
        scale_pga = np.zeros(R.shape)
        ycor_surface = np.zeros(R.shape)
        ampli = np.zeros(R.shape)
        if mode == 'median':
            y = pow(10., (a1 + a2 * M + a3 * M ** 2 + a4 * M ** 3 + a5 * M ** 4 + a6 * M ** 5 + a7 * M ** 6
                          + (a8 + a9 * M + a10 * M * M + a11 * M * M * M) * d
                          + (a12 + a13 * M + a14 * M * M + a15 * M * M * M) * (d ** 2)
                          + (a16 + a17 * M + a18 * M * M + a19 * M * M * M) * (d ** 3)
                          + (a20 + a21 * M + a22 * M * M + a23 * M * M * M) * (d ** 4)))
        elif mode == '84':
            y = pow(10., (a1 + a2 * M + a3 * M ** 2 + a4 * M ** 3 + a5 * M ** 4 + a6 * M ** 5 + a7 * M ** 6
                          + (a8 + a9 * M + a10 * M * M + a11 * M * M * M) * d
                          + (a12 + a13 * M + a14 * M * M + a15 * M * M * M) * (d ** 2)
                          + (a16 + a17 * M + a18 * M * M + a19 * M * M * M) * (d ** 3)
                          + (a20 + a21 * M + a22 * M * M + a23 * M * M * M) * (d ** 4) + 0.2017))
        if sn:
            # Correction to local basement at Beznau
            ycor = y * 0.83771575
            # Surface amplification as a function of M and log10pga in units of g
            # simplified fitting swissnuclear data ...
            scale_pga = np.log10(ycor / 981)
            if M <= 5.25:
                ampli = -0.57 * scale_pga ** 2 - 1.7 * scale_pga + 1.2
                ampli = np.where(ycor / 981 < 0.05, 2.4469, ampli)
            elif M <= 6.25:
                ampli = -0.63 * scale_pga ** 2 - 1.6 * scale_pga + 1.2
                ampli = np.where(ycor / 981 < 0.05, 2.2153, ampli)
            else:
                ampli = -0.64 * scale_pga ** 2 - 1.7 * scale_pga + 1.2;
                ampli = np.where(ycor / 981 < 0.05, 2.3284, ampli)
            ycor_surface = ycor * ampli
        else:
            ycor_surface = y
        return R, ycor_surface

if __name__ == '__main__':
    sg = SwissGmpe()
    distance = np.arange(30, 101, 1)
    for mag in np.arange(3.25, 6.75, 0.5):
        R, ycor_surface = sg.get_pga(mag, distance, 'foreland', 'median')
        plt.loglog(R, ycor_surface)
        R, ycor_surface = sg.get_pga(mag, np.array([distance[0]]), 'foreland', 'median')
        plt.loglog(R, ycor_surface, 'ko')
        R, ycor_surface = sg.get_pga(mag, np.array([distance[-1]]), 'foreland', 'median')
        plt.loglog(R, ycor_surface, 'ko')

    # plt.legend()
    plt.xlim(30, 100)
    plt.ylim(0.1, 100)
    plt.show()




