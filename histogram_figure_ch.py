#!/usr/bin/env python
"""
Combine all delay distributions in one figure.
Created on Nov 26, 2013

@author: behry
"""
import os
import numpy as np
import matplotlib
#matplotlib.use('WXAgg')
try: # added by FM
    matplotlib.use('WXAgg')
    import wx
except:
    print "WX package for Python not installed"
import matplotlib.pyplot as plt
import json
from matplotlib import rcParams
from timeline import TimeLine
from scipy.stats import scoreatpercentile

rcParams['figure.subplot.left'] = 0.03
rcParams['figure.subplot.right'] = 0.97
rcParams['figure.subplot.top'] = 1.00
rcParams['figure.subplot.bottom'] = 0.05
rcParams['figure.subplot.hspace'] = 0.00
rcParams['figure.subplot.wspace'] = 0.00
rcParams['legend.fontsize'] = 12
rcParams['xtick.major.size'] = 6


class HistPlot:

    def __init__(self, data, ax):
        self.data = data
        self.ax = ax
        self.bins = np.arange(0, 20, 0.5)
        self.colors = ['blue']
        self.normed = True
        self.label = []

    def text(self):
        starty = 0.5
        startx = 0.66
        dy = 0.08
        alldel = []
        meds = []
        datamin = self.bins[0]
        datamax = self.bins[-1]
        if len(self.data) < 2:
            self.colors = ['black']
        for _i, _d in enumerate(self.data):
            alldel += _d
            med = np.ma.median(np.ma.masked_where((_d < datamin) | (_d > datamax), _d))
            # med = np.median(_d)
            meds.append(med)
            self.ax.text(startx, starty, 'Median: %.1f s' % (med),
                    horizontalalignment='left',
                    transform=self.ax.transAxes, color=self.colors[_i])
            starty -= dy
        meds.append(np.median(alldel))
        percentile16 = scoreatpercentile(alldel, 16)
        percentile84 = scoreatpercentile(alldel, 84)
        if len(self.data) > 1:
            self.ax.text(startx, starty, 'Median combined: %.1f s' % (meds[-1]),
                         horizontalalignment='left',
                         transform=self.ax.transAxes, color='black')
            starty -= dy
        self.ax.text(startx, starty, '16th percentile: %.1f s' % (percentile16),
                horizontalalignment='left',
                transform=self.ax.transAxes, color='black')
        starty -= dy
        self.ax.text(startx, starty, '84th percentile: %.1f s' % (percentile84),
                horizontalalignment='left',
                transform=self.ax.transAxes, color='black')

    def annotate(self, lbl):
        startx = 0.12
        starty = 0.9
        self.ax.text(startx, starty, lbl, horizontalalignment='left',
                     transform=self.ax.transAxes, color='black', fontsize=14,
                     weight='bold')

    def delay_description(self, desc):
        startx = 0.35
        starty = 0.72
        self.ax.text(startx, starty, desc, horizontalalignment='left',
                     transform=self.ax.transAxes, color='black',
                     fontsize=13)

    def hist(self):
        self.ax.hist(self.data, bins=self.bins, color=self.colors,
                     rwidth=1.0, normed=self.normed, histtype='barstacked',
                     label=self.label)
        self.ax.set_xticks([0, 5, 10, 15, 20])
        self.ax.set_xticklabels([])
        self.ax.set_yticks([])
        self.ax.legend()

    def plot(self, lbl=''):
        self.hist()
        self.text()
        self.annotate(lbl)



def data_latencies(ax, fn):
    # Data latency
    a = np.load(fn)
    datlat = a['stations']
    chdatlat = []
    frgdatlat = []
    for _e in datlat:
        stnm = _e[0]['name']
        val = _e[0]['median']
        net, stat = stnm.split('.')
        if net == 'CH' or net == 'S' or net == '8D':
            chdatlat.append(val)
        else:
            frgdatlat.append(val)

    hp = HistPlot([chdatlat, frgdatlat], ax)
    hp.colors = ['blue', 'green']
    hp.label = ['Swiss stations', 'Foreign stations']
    hp.plot('a')
    hp.delay_description(r'Data latency (' + r'$\Delta t_{l,i}$' + r')')
    ax.set_title('Switzerland', fontsize=16)

def envelope_delays(ax, fn, perstation=False):
    fh = open(fn)
    envdel = json.load(fh)
    fh.close()
    chdelayse = []
    frgdelayse = []
    for _k in envdel.keys():
        net, stat = _k.split('.')
        if net == 'CH' or net == 'S' or net == '8D':
            if perstation:
                chdelayse.append(np.median(envdel[_k]))
            else:
                chdelayse += envdel[_k]
        else:
            if perstation:
                frgdelayse.append(np.median(envdel[_k]))
            else:
                frgdelayse += envdel[_k]

    hp = HistPlot([chdelayse, frgdelayse], ax)
    hp.colors = ['blue', 'green']
    hp.label = ['Swiss stations', 'Foreign stations']
    hp.plot('c')
    hp.delay_description(r'Waveform delay (' + r'$\Delta t_{w,i}$' + r')')



def pick_delays(ax, fn, perstation=False):
    fh = open(fn)
    pkdel = json.load(fh)
    fh.close()
    chdelays = []
    frgdelays = []
    for _k in pkdel.keys():
        net, stat = _k.split('.')
        if net == 'CH' or net == 'S' or net == '8D':
            if perstation:
                chdelays.append(np.median(pkdel[_k]))
            else:
                chdelays += pkdel[_k]
        else:
            if perstation:
                frgdelays.append(np.median(pkdel[_k]))
            else:
                frgdelays += pkdel[_k]

    hp = HistPlot([chdelays, frgdelays], ax)
    hp.colors = ['blue', 'green']
    hp.label = ['Swiss stations', 'Foreign stations']
    hp.plot('e')
    hp.delay_description(r'Trigger delay (' + r'$\Delta t_{t,i}$' + r')')


def associator_delays(ax, fn):
    a = np.load(fn)
    octdel = a['delays']
    hp = HistPlot([octdel.tolist()], ax)
    hp.plot('g')
    hp.delay_description(r'Associator delay (' + r'$\Delta t_a$' + r')')


def magnitude_delays(ax, fn):
    a = np.load(fn)
    magdel = a['delays']
    hp = HistPlot([magdel.tolist()], ax)
    hp.plot('i')
    hp.delay_description(r'Magnitude delay (' + r'$\Delta t_m$' + r')')


def histogram_ch(fig, sp, dsp):
    rows, cols, idx = sp
    # get all the delay distributions
    datadir = './data'
    datlatfn = os.path.join(datadir, 'data_delays_ch.npz')
    pkdelfn = os.path.join(datadir, 'single_station_pk_delays_ch.txt')
    envdelfn = os.path.join(datadir, 'single_station_env_delays_ch.txt')
    octdelfn = os.path.join(datadir, 'origin_estimation_delay_ch.npz')
    magdelfn = os.path.join(datadir, 'magnitude_estimation_delay.npz')

    perstation = True
    # Data latencies
    ax = fig.add_subplot(rows, cols, idx)
    data_latencies(ax, datlatfn)

    # Envelope delays
    idx += dsp
    ax = fig.add_subplot(rows, cols, idx)
    envelope_delays(ax, envdelfn, perstation=perstation)

    # Pick delays
    idx += dsp
    ax = fig.add_subplot(rows, cols, idx)
    pick_delays(ax, pkdelfn, perstation=perstation)

    # Associator delays
    idx += dsp
    ax = fig.add_subplot(rows, cols, idx)
    associator_delays(ax, octdelfn)

    # Magnitude delays
    idx += dsp
    ax = fig.add_subplot(rows, cols, idx)
    magnitude_delays(ax, magdelfn)
    ax.set_xticklabels(['%d' % i for i in ax.get_xticks()])
    ax.set_xlabel('Time since origin time [s]', fontsize=14)

if __name__ == '__main__':
    fout = './plots/delay_histograms_ch.png'
    fig = plt.figure(figsize=(8.267, 11.692))
    histogram_ch(fig, (5, 1, 1), 1)
    plt.savefig(fout, dpi=300,
                bbox_inches='tight',
                transparent=True)
    plt.show()



