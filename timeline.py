#!/usr/bin/env python
"""
Created on Nov 26, 2013

@author: behry
"""
import sys
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
# rcParams['savefig.format'] = 'ps'
rcParams['figure.subplot.bottom'] = 0.42
#     rcParams['axes.edgecolor'] = 'None'
#     rcParams['axes.facecolor'] = 'None'
#     rcParams['axes.axisbelow'] = 'True'
#     rcParams['xtick.direction'] = 'out'


class TimeLine:

    def __init__(self, ax, fontsize=10, text=True, hlt=''):
        self.cls = [(138 / 255., 138 / 255., 230 / 255.),
                    (181 / 255., 204 / 255., 112 / 255.),
                    (145 / 255., 243 / 255., 235 / 255.),
                    (172 / 255., 243 / 255., 133 / 255.),
                    (217 / 255., 163 / 255., 163 / 255.),
                    (217 / 255., 163 / 255., 220 / 255.)]
        self.ax = ax
        self.height = 0.5
        self.bx_offst = 0.3
        self.fontsize = fontsize
        self.xmin = 0
        self.xmax = 15.4
        self.ymin = 0
        self.ymax = 0
        self.torg = 0
        self.text = text
        self.hlt = hlt
        self.counter = 'i'

    def add_preproc(self, xstart=0.0, ystart=0.0, tp=7.5, tl=1.5, tt=1., tw=3.,
                    top=False, bottom=False, label='', version='detailed'):
        self.ymax = max(ystart + 3 * self.height, self.ymax)
        self.xmin = min(xstart, self.xmin)
        self.torg = max(tp + tl + tt, self.torg)
        txt_offst = 0.15
        self.ax.bar(xstart, self.height, tp, ystart + self.height, orientation='horizontal',
               color=self.cls[0])
        if version == 'detailed':
            p = self.ax.bar(tp, self.height, tl, ystart + 2 * self.height,
                            orientation='horizontal', color='gray')
            self.ax.bar(tp + tl, self.height, tw, ystart + 2 * self.height,
                        orientation='horizontal', color=self.cls[1])
            self.ax.bar(tp, self.height, tl, ystart, orientation='horizontal', color='gray')
            self.ax.bar(tp + tl, self.height, tt, ystart, orientation='horizontal',
                        color=self.cls[2])
        elif version == 'simple':
            p = self.ax.bar(tp, self.height, tl + tw, ystart + 2 * self.height,
                            orientation='horizontal', color=self.cls[1])
            self.ax.bar(tp, self.height, tl + tt, ystart,
                        orientation='horizontal', color=self.cls[2])

        if self.text:
            self.ax.text((tp - xstart) / 2., ystart + self.height + txt_offst,
                         r'$\mathrm{\Delta t_{p,%s}}$' % self.counter,
                         fontsize=self.fontsize,
                         fontweight='bold', horizontalalignment='center')
            if version == 'detailed':
                self.ax.text(xstart + tp + tl / 2., ystart + 2 * self.height + txt_offst,
                        r'$\mathrm{\Delta t_{l,%s}}$' % self.counter,
                        fontsize=self.fontsize,
                        fontweight='bold', horizontalalignment='center')
                self.ax.text(xstart + tp + tl + tw / 2.,
                        ystart + 2 * self.height + txt_offst,
                        r'$\mathrm{\Delta t_{pw} + \Delta t_{win}}$',
                        fontsize=self.fontsize, fontweight='bold',
                        horizontalalignment='center')
                self.ax.text(xstart + tp + tl / 2., ystart + txt_offst,
                        r'$\mathrm{\Delta t_{l,%s}}$' % self.counter,
                        fontsize=self.fontsize,
                        fontweight='bold', horizontalalignment='center')
                self.ax.text(xstart + tp + tl + tt / 2., ystart + txt_offst,
                        r'$\mathrm{\Delta t_{pk}}$', fontsize=self.fontsize, fontweight='bold',
                        horizontalalignment='center')

                self.ax.text(xstart - txt_offst , ystart + self.height + txt_offst,
                             label, fontsize=self.fontsize, fontweight='bold',
                             horizontalalignment='right')
                offset = 1.
                if top:
                    self.ax.annotate(r'$\mathrm{\Delta t_{w,%s}}$' % self.counter,
                                     (tp + (tl + tw) / 2., ystart + 3 * self.height + 0.1),
                                     (tp + (tl + tw) / 2. , ystart + 5 * self.height),
                                     ha="center", va="center",
                                     fontsize=self.fontsize,
                                     arrowprops=dict(arrowstyle='-[', shrinkA=1,
                                                     shrinkB=20, fc='w', ec='k',
                                                     mutation_scale=105))
                if bottom:
                    self.ax.annotate(r'$\mathrm{\Delta t_{t,%s}}$' % self.counter,
                                     (tp + (tl + tt) / 2., 0),
                                     (tp + (tl + tt) / 2., ystart - 2 * self.height),
                                     ha="center", va="center",
                                     fontsize=self.fontsize,
                                     arrowprops=dict(arrowstyle='-[', shrinkA=5,
                                                     shrinkB=15, fc='w', ec='k',
                                                     mutation_scale=52))
            elif version == 'simple':
                self.ax.text(xstart - txt_offst , ystart + self.height + txt_offst,
                             label, fontsize=self.fontsize, fontweight='bold',
                             horizontalalignment='right')
                self.ax.text(xstart + tp + (tl + tw) / 2., ystart + 2 * self.height + txt_offst,
                             r'$\mathrm{\Delta t_{w,%s}}$' % self.counter,
                             fontsize=self.fontsize,
                             fontweight='bold', horizontalalignment='center')
                self.ax.text(xstart + tp + (tl + tt) / 2., ystart + txt_offst,
                             r'$\mathrm{\Delta t_{t,%s}}$' % self.counter,
                             fontsize=self.fontsize, fontweight='bold',
                             horizontalalignment='center')

            if self.hlt == 'tp':
                self.ax.bar(xstart - self.bx_offst,
                            self.height + 2 * self.bx_offst,
                            tp + 2 * self.bx_offst,
                            ystart + self.height - self.bx_offst,
                            edgecolor='black', facecolor='None', linewidth=2.0)
            elif self.hlt == 'tl':
                self.ax.bar(xstart + tp - self.bx_offst,
                            3 * self.height + 2 * self.bx_offst,
                            tl + 2 * self.bx_offst,
                            ystart - self.bx_offst, edgecolor='black',
                            facecolor='None', linewidth=2.0)
            elif self.hlt == 'tw':
                self.ax.bar(xstart + tp - self.bx_offst,
                            self.height + 2 * self.bx_offst,
                            tl + tw + 2 * self.bx_offst,
                            ystart + 2 * self.height - self.bx_offst,
                            edgecolor='black',
                            facecolor='None', linewidth=2.0)
            elif self.hlt == 'tt':
                self.ax.bar(xstart + tp - self.bx_offst,
                            self.height + 2 * self.bx_offst,
                            tl + tt + 2 * self.bx_offst,
                            ystart - self.bx_offst, edgecolor='black',
                            facecolor='None', linewidth=2.0)

    def add_post_proc(self, ta, tm, td, txt_offst=0.15):
        self.xmax = max(self.torg + ta + tm + td + self.bx_offst, self.xmax)
        self.ax.vlines(self.torg, self.ymin - self.height,
                       self.ymax + self.height,
                       linestyles='dashed')
        y0 = (self.ymax - self.ymin) / 2. - self.height / 2.
        self.ax.bar(self.torg, self.height, ta, y0,
                    orientation='horizontal', color=self.cls[3])
        self.ax.bar(self.torg + ta, self.height, tm, y0,
               orientation='horizontal', color=self.cls[4])
        self.ax.bar(self.torg + ta + tm, self.height, td, y0,
               orientation='horizontal', color=self.cls[5])
        if self.text:
            self.ax.text(self.torg + ta / 2., y0 + txt_offst,
                         r'$\mathrm{\Delta t_a}$', fontsize=self.fontsize,
                         fontweight='bold', horizontalalignment='center')
            self.ax.text(self.torg + ta + tm / 2.,
                         y0 + txt_offst, r'$\mathrm{\Delta t_m}$',
                         fontsize=self.fontsize, fontweight='bold',
                         horizontalalignment='center')
            self.ax.text(self.torg + ta + tm + td / 2.,
                         y0 + txt_offst, r'$\mathrm{\Delta t_d}$',
                         fontsize=self.fontsize, fontweight='bold',
                         horizontalalignment='center')
            self.ax.annotate(r'$\mathrm{\Delta t_{origin}}$',
                             (self.torg + ta, y0 + 1 * self.height),
                             (self.torg + ta , y0 + 2 * self.height),
                             ha="center", va="center",
                             fontsize=self.fontsize,
                             arrowprops=dict(arrowstyle='-|>', shrinkA=1,
                                             shrinkB=1, fc='w', ec='k'))
            self.ax.annotate(r'$\mathrm{\Delta t_{alert}}$',
                        (self.torg + ta + tm + td, y0 + 1 * self.height),
                        (self.torg + ta + tm + td, y0 + 2 * self.height),
                        ha="center", va="center",
                        fontsize=self.fontsize,
                        arrowprops=dict(arrowstyle='-|>', shrinkA=1, shrinkB=1,
                                        fc='w', ec='k'))
        if self.hlt == 'ta':
            self.ax.bar(self.torg - self.bx_offst,
                        self.height + 2 * self.bx_offst,
                        ta + 2 * self.bx_offst,
                        y0 - self.bx_offst, edgecolor='black',
                        facecolor='None', linewidth=2.0)
        elif self.hlt == 'tm':
            self.ax.bar(self.torg + ta - self.bx_offst,
                        self.height + 2 * self.bx_offst,
                        tm + 2 * self.bx_offst,
                        y0 - self.bx_offst,
                        edgecolor='black', facecolor='None', linewidth=2.0)

    def add_timeaxis(self):
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.hlines(-2 * self.height - 2 * self.bx_offst, 0, 15)
        self.ax.set_xlabel('Time since origin time [s]', fontsize=self.fontsize)

    def set_limits(self):
        self.ax.set_ylim(-2 * self.height - 2 * self.bx_offst,
                         self.ymax + 2 * self.height + 2 * self.bx_offst)
        self.ax.set_xlim(-2 * self.bx_offst, max(self.xmax, 15.4))
        self.ax.set_yticks([])
        self.ax.set_axis_bgcolor('None')
        self.ax.set_frame_on(False)

    def vdots(self):
        self.ax.vlines(self.xmin - 2. * self.bx_offst,
                       self.ymin + 2. * self.height,
                       self.ymin + 5. * self.height, colors='k',
                       linestyles='dotted', linewidth=5)


if __name__ == '__main__':
    if True:
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111)
        tl = TimeLine(ax, fontsize=12)
        tl.counter = '1'
        tl.add_preproc(tp=3.5, tl=2.0, ystart=4.0, top=True, label='Station 1',
                       version='simple')
        tl.counter = '2'
        tl.add_preproc(tp=5.5, tl=1.0, ystart=2.0, label='Station 2',
                       version='simple')
        tl.counter = 'n'
        tl.add_preproc(bottom=True, label='Station n', version='simple')
        tl.add_post_proc(1.2, 1., 0.8, txt_offst=0.15)
        tl.add_timeaxis()
        tl.vdots()
        tl.set_limits()
        plt.savefig('/home/behry/uni/eew/SSA_2015/pix/timeline.png', dpi=300,
                    bbox_inches='tight')
        plt.show()
    if False:
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        tl = TimeLine(ax, text=False, hlt='tm')
        tl.add_preproc()
        tl.torg = 12.5
        tl.add_post_proc(1.2, 1., 0.8, txt_offst=0.15)
        tl.set_limits()
        ax.set_axis_off()
        plt.show()

