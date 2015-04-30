#!/usr/bin/env python
"""
Compute whether an event is located inside or outside of Switzerland.
Created on Sep 10, 2012

@author: behry
"""


from mpl_toolkits.basemap import Basemap
import shapefile
import pyproj
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, MultiPoint
from shapely.prepared import prep

class EventCh:
    def __init__(self):
        self.sf = shapefile.Reader('/home/behry/uni/data/shapefiles/shp/gg25_a_v2_Dissolve')
        shapes = self.sf.shapes()
        self.points = np.array(shapes[0].points)
        # Convert the Swiss Coordinate System into lat/lon
        self.p = pyproj.Proj(proj="somerc", lat_0='46d57\'8.660"N', lon_0='7d26\'22.500"E',
                        ellps='bessel', x_0=600000, y_0=200000, k_0=1.)
        self.polygon = Polygon(self.points)
        self.prepared_polygon = prep(self.polygon)
        self.lat1 = 45.5
        self.lat2 = 48
        self.lon1 = 5
        self.lon2 = 12

    def point_in_polygon(self, lon, lat):
        p1x, p1y = self.p(lon, lat)
        return np.array(map(self.prepared_polygon.contains, map(Point, zip(p1x, p1y))))

    def plot_polygon(self, point=None):
        m = Basemap(projection='merc', llcrnrlat=self.lat1,
                    urcrnrlat=self.lat2, llcrnrlon=self.lon1,
                    urcrnrlon=self.lon2, lat_ts=47, resolution='i')
        lon, lat = self.p(self.points[:, 0], self.points[:, 1], inverse=True)
        x, y = m(lon, lat)
        m.drawcountries(linewidth='1.0')
        m.drawstates()
        m.drawcoastlines()
        m.plot(x, y, 'r--')
        if point is not None:
            xp, yp = m(point[0], point[1])
            m.plot(xp, yp, 'bo')
        plt.show()

    def plot_polygon_xy(self, point=None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        patch1 = PolygonPatch(self.polygon, facecolor='#6699cc',
                              edgecolor='#6699cc', alpha=0.5, zorder=2)
        ax.add_patch(patch1)
        ax.plot(self.points[:, 0], self.points[:, 1], 'r--')
        if point is not None:
            p1x, p1y = self.p(point[0], point[1])
            plt.plot(p1x, p1y, 'bo')
        plt.show()

class EventCA(EventCh):

    def __init__(self):
        # filterfn = '/home/behry/workspace/eew/reports/data/california_filter_revised.txt'
        filterfn = '/home/behry/workspace/eew/reports/data/cafilter_dm_review.txt'
        lon, lat = np.loadtxt(filterfn, unpack=True)
        self.p = pyproj.Proj(proj="aea", lat_0='0.00', lon_0='-120',
                             lat_1='34.00', lat_2='40.50',
                             ellps='GRS80', x_0=0.000, y_0=-4000000.000,
                             units='m', datum='NAD83')
        px, py = self.p(lon, lat)
        self.points = np.vstack((px, py)).T
        self.polygon = Polygon(zip(px, py))
        self.prepared_polygon = prep(self.polygon)
        self.lat1 = 31
        self.lat2 = 43
        self.lon1 = -126
        self.lon2 = -114

class EventSoCal(EventCh):

    def __init__(self):
        # filterfn = '/home/behry/workspace/eew/reports/data/california_filter_revised.txt'
        filterfn = '/home/behry/workspace/eew/reports/data/SoCalfilter.txt'
        lon, lat = np.loadtxt(filterfn, unpack=True)
        self.p = pyproj.Proj(proj="aea", lat_0='0.00', lon_0='-120',
                             lat_1='34.00', lat_2='40.50',
                             ellps='GRS80', x_0=0.000, y_0=-4000000.000,
                             units='m', datum='NAD83')
        px, py = self.p(lon, lat)
        self.points = np.vstack((px, py)).T
        self.polygon = Polygon(zip(px, py))
        self.prepared_polygon = prep(self.polygon)
        self.lat1 = 31
        self.lat2 = 43
        self.lon1 = -126
        self.lon2 = -114

class EventCHP(EventCh):

    def __init__(self):
        filterfn = '/home/behry/workspace/eew/reports/data/swiss_polygon.txt'
        lon, lat = np.loadtxt(filterfn, unpack=True)
        self.p = pyproj.Proj(proj="somerc", lat_0='46d57\'8.660"N',
                             lon_0='7d26\'22.500"E',
                             ellps='bessel', x_0=600000, y_0=200000, k_0=1.)
        px, py = self.p(lon, lat)
        self.points = np.vstack((px, py)).T
        self.polygon = Polygon(zip(px, py))
        self.prepared_polygon = prep(self.polygon)
        self.lat1 = 45.5
        self.lat2 = 48
        self.lon1 = 5
        self.lon2 = 12

class EventBrawley(EventCh):

    def __init__(self):
        filterfn = '/home/behry/workspace/eew/reports/data/brawley_polygon.txt'
        lon, lat = np.loadtxt(filterfn, unpack=True)
        self.p = pyproj.Proj(proj="aea", lat_0='0.00', lon_0='-120',
                             lat_1='34.00', lat_2='40.50',
                             ellps='GRS80', x_0=0.000, y_0=-4000000.000,
                             units='m', datum='NAD83')
        px, py = self.p(lon, lat)
        self.points = np.vstack((px, py)).T
        self.polygon = Polygon(zip(px, py))
        self.prepared_polygon = prep(self.polygon)
        self.lat1 = 31
        self.lat2 = 43
        self.lon1 = -126
        self.lon2 = -114


if __name__ == '__main__':
    if False:
        t = EventCh()
        # Zurich
        point1 = np.array([8.540, 47.372])
        # Waldshut-Tiengen
        point2 = np.array([8.272705, 47.642261])

        print t.point_in_polygon([point1[0]], [point1[1]])
        print t.point_in_polygon([point2[0]], [point2[1]])
        t.plot_polygon()
        t.plot_polygon_xy()

    if True:
        t = EventSoCal()
        t.plot_polygon()

    if False:
        t = EventCA()
        # Reno
        point1 = np.array([-119.793755, 39.533703])
        # Truckee
        point2 = np.array([-120.187889, 39.330049])
        # Eureka
        point3 = np.array([-124.166897, 40.794058])
        # multiple points
        lats = [39.533703, 39.330049, 40.794058]
        lons = [-119.793755, -120.187889, -124.166897]
        print t.point_in_polygon([point1[0]], [point1[1]])
        print t.point_in_polygon([point2[0]], [point2[1]])
        print t.point_in_polygon([point3[0]], [point3[1]])
        print t.point_in_polygon(lons, lats)
        t.plot_polygon_xy(point3)
        t.plot_polygon(point3)

    if False:
        t = EventCHP()
        t.plot_polygon()

    if False:
        t = EventBrawley()
        t.plot_polygon()
