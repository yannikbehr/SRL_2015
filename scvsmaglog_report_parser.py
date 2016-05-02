#!/usr/bin/env python
"""
Parse reports generated by SeisComp3's scvsmaglog module.
Created on Mar 10, 2013

@author: behry
"""
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth # obspy.core.util.gps2DistAzimuth soon to be deprecated
import numpy as np

def parser(fn, trueot=None, trueloc=None, format='old'):
    fh = open(fn)
    contents = []
    # skip header lines
    fh.readline()
    fh.readline()
    for l in fh.readlines():
        l = l.rstrip()
        a = l.split('|')
        mag = float(a[0])
        lat = float(a[1])
        lon = float(a[2])
        if format == 'old':
            dep = float(a[3])
            ct = UTCDateTime(a[4])
            ot = UTCDateTime(a[5])
            likeh = float(a[-1])
        else:
            dep = float(a[4])
            ct = UTCDateTime(a[5])
            ot = UTCDateTime(a[6])
            likeh = float(a[7])
        if trueot is None:
            tdiff = float(ct - ot)
        else:
            tdiff = float(ct - trueot)
        if trueloc is None:
            ddiff = 0
        else:
            dist, az, baz = gps2DistAzimuth(lat, lon, *trueloc)
            ddiff = dist / 1000.
        if format != 'old':
            nstorig = int(a[-2])
            nstmag = int(a[-1])
            temp = [mag, lat, lon, dep, ct, ot, tdiff, ddiff,
                    likeh, nstorig, nstmag]
        else:
            temp = [mag, lat, lon, dep, ct, ot, tdiff, ddiff, likeh]
        contents.append(temp)

    return np.array(contents)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    fn = '/home/behry/workspace/eew/reports/data/scvsmag_report_zug.txt'
    cts = parser(fn)
    plt.plot(cts[:, 6], cts[:, 0])
    plt.show()
