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
from histogram_figure_ch import HistPlot

rcParams['figure.subplot.left'] = 0.03
rcParams['figure.subplot.right'] = 0.97
rcParams['figure.subplot.top'] = 1.00
rcParams['figure.subplot.bottom'] = 0.05
rcParams['figure.subplot.hspace'] = 0.00
rcParams['figure.subplot.wspace'] = 0.00
rcParams['legend.fontsize'] = 12


class HistPlotCA(HistPlot):

    def __init__(self, data, ax, fig):
        HistPlot.__init__(self, data, ax)
        self.fig = fig


    def inset(self, hlt):
        Bbox = matplotlib.transforms.Bbox.from_bounds(.22, .55, .48, .35)
        trans = self.ax.transAxes + self.fig.transFigure.inverted()
        l, b, w, h = matplotlib.transforms.TransformedBbox(Bbox, trans).bounds
        axins = self.fig.add_axes([l, b, w, h])
        tl = TimeLine(axins, text=False, hlt=hlt)
        tl.add_preproc()
        tl.torg = 12.5
        tl.add_post_proc(1.2, 1., 0.8, txt_offst=0.15)
        tl.set_limits()
        axins.set_axis_off()


def load_delays(fn, perstation=False):
    fh = open(fn)
    delays = json.load(fh)
    fh.close()
    bkdelays = []
    mpdelays = []
    cadelays = []
    for _k in delays.keys():
        net, stat = _k.split('.')
        if net == 'NP' or net == 'NC':
            if perstation:
                mpdelays.append(np.median(delays[_k]))
            else:
                mpdelays += delays[_k]
        elif net == 'CI' or net == 'AZ':
            if perstation:
                cadelays.append(np.median(delays[_k]))
            else:
                cadelays += delays[_k]
        elif net == 'BK':
            if perstation:
                bkdelays.append(np.median(delays[_k]))
            else:
                bkdelays += delays[_k]
        else:
            print "Unknown network %s" % net
    return bkdelays, mpdelays, cadelays


def data_latencies(fig, ax, fn, perstation=False):
    bkdelays, mpdelays, cadelays = load_delays(fn, perstation=perstation)
    hp = HistPlotCA([cadelays, bkdelays, mpdelays], ax, fig)
    hp.colors = ['blue', 'green', 'red']
    hp.label = ['Caltech', 'Berkeley', 'Menlo Park']
    hp.plot('b')
    hp.inset('tl')
    ax.set_title('California', fontsize=16)


def envelope_delays(fig, ax, fn, perstation=False):
    bkdelays, mpdelays, cadelays = load_delays(fn, perstation=perstation)
    hp = HistPlotCA([cadelays, bkdelays, mpdelays], ax, fig)
    hp.colors = ['blue', 'green', 'red']
    hp.label = ['Caltech', 'Berkeley', 'Menlo Park']
    hp.plot('d')
    hp.inset('tw')


def pick_delays(fig, ax, fn, perstation=False):
    bkdelays, mpdelays, cadelays = load_delays(fn, perstation=perstation)
    hp = HistPlotCA([cadelays, bkdelays, mpdelays], ax, fig)
    hp.colors = ['blue', 'green', 'red']
    hp.label = ['Caltech', 'Berkeley', 'Menlo Park']
    hp.plot('f')
    hp.inset('tt')


def associator_delays(fig, ax, fn):
    a = np.load(fn)
    octdel = a['acdl']
    hp = HistPlotCA([octdel.tolist()], ax, fig)
    hp.plot('h')
    hp.inset('ta')


def magnitude_delays(fig, ax, fn):
    a = np.load(fn)
    magdel = a['mct']
    hp = HistPlotCA([magdel.tolist()], ax, fig)
    hp.plot('j')
    hp.inset('tm')


def histogram_ca(fig, sp, dsp):
    rows, cols, idx = sp
    # get all the delay distributions
    datadir = './data'
    datlatfn = os.path.join(datadir, 'data_delays_ca.txt')
    pkdelfn = os.path.join(datadir, 'single_station_pk_delays_ca.txt')
    envdelfn = os.path.join(datadir, 'single_station_env_delays_ca.txt')
    octdelfn = os.path.join(datadir, 'eewvs_origin_delays.npz')
    magdelfn = os.path.join(datadir, 'eewvs_magnitude_delays.npz')

    perstation = True
    # Data latencies
    ax = fig.add_subplot(rows, cols, idx)
    data_latencies(fig, ax, datlatfn, perstation=perstation)

    # Envelope delays
    idx += dsp
    ax = fig.add_subplot(rows, cols, idx)
    envelope_delays(fig, ax, envdelfn, perstation=perstation)

    # Pick delays
    idx += dsp
    ax = fig.add_subplot(rows, cols, idx)
    pick_delays(fig, ax, pkdelfn, perstation=perstation)

    # Associator delays
    idx += dsp
    ax = fig.add_subplot(rows, cols, idx)
    associator_delays(fig, ax, octdelfn)

    # Magnitude delays
    idx += dsp
    ax = fig.add_subplot(rows, cols, idx)
    magnitude_delays(fig, ax, magdelfn)

    ax.set_xticklabels(['%d' % i for i in ax.get_xticks()])
    ax.set_xlabel('Time since origin time [s]', fontsize=14)


if __name__ == '__main__':
    fout = './plots/delay_histograms_ca.png'
    fig = plt.figure(figsize=(8.267, 11.692))
    histogram_ca(fig, (5, 1, 1), 1)
    plt.savefig(fout, dpi=300,
                bbox_inches='tight',
                transparent=True)
    plt.show()
