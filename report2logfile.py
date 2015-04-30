#!/usr/bin/env python
"""
Find the log file that corresponds to a report and extract more details
on the particular event from the log file.
Created on Oct 31, 2013

@author: behry
"""

from obspy import UTCDateTime
import os
import sys
import subprocess as sp
from eewvs_log_parser import read_eewvs_log
from reports_parser import ReportsParser
import re
import json
import pyproj


def scp_copy_logfile(filein, fileout, datdir):
    if not os.path.isfile(fileout):
        sshcommand = "scp -i /home/behry/.ssh/id_rsa_havana "
        sshcommand += "eewvs@havana.ethz.ch:%s %s.bz2" % (filein, fileout)
        print sshcommand
        cwd = os.getcwd()
        sp.call(sshcommand, shell=True)
        os.chdir(datdir)
        sp.call("bunzip2 %s.bz2" % fileout, shell=True)
        os.chdir(cwd)
        if not os.path.isfile(fileout):
            return False
    return True


def event_in_southern_california(publicID):
    if publicID.startswith('CI'):
        return True
    return False


def get_event_info(dt, reportfile, datdir, publicID):
    print reportfile, publicID
    try:
        rid = int(re.match(r'eewvs-report-\d+-\d+-(\d+)\.xml', reportfile).group(1))
    except AttributeError:
        print reportfile
        raise AttributeError()

    logfile = "eewvs-%s.log" % dt.strftime("%Y%m%d")
    logfile1 = '/home/projects/eewvs/wwwhome/cinnabar/reports/'
    logfile1 += '%d/%02d/%s.bz2' % (dt.year, dt.month, logfile)
    localfile1 = "%s/%s.cinnabar" % (datdir, logfile)
    logfile2 = '/home/projects/eewvs/wwwhome/berkeley/reports/'
    logfile2 += '%d/%02d/%s.bz2' % (dt.year, dt.month, logfile)
    localfile2 = "%s/%s.berkeley" % (datdir, logfile)
    logfile3 = '/home/projects/eewvs/wwwhome/menlopark/reports/'
    logfile3 += '%d/%02d/%s.bz2' % (dt.year, dt.month, logfile)
    localfile3 = "%s/%s.menlopark" % (datdir, logfile)

    evid = None
    for _f1, _f2 in [(logfile3, localfile3),
                     (logfile2, localfile2),
                     (logfile1, localfile1)]:
        if not scp_copy_logfile(_f1, _f2, datdir):
            continue
        ev, pk, dummy = read_eewvs_log(_f2, format='old')
        for _id in ev.keys():
            try:
                if ev[_id]['report'] == reportfile:
                    evid = _id
                    picks = pk
            except:
                continue
        if evid is not None:
            break
    if evid is None:
        print "Corresponding log-file not found for event %s" % (publicID)
        return []

    # Consistency check
    if not evid == rid:
        raise Exception("Event IDs don't match! (evid=%d, rid=%d" % (evid, rid))
    try:
        estno, nnodes, mag, evlat, evlon, depth, rtime, \
        vtime, ts, mguess, merr, merrthresh, pickerr = ev[evid]['estimates'][0]
    except:
        raise Exception('No valid estimates found for %s: %d' % \
                        (reportfile, ev[evid]['cnt']))

    # Find the nodes update just before the first estimate
    _vts = UTCDateTime(vtime)
    _nodes = 1
    inst = '%04d' % (_nodes)
    if _nodes > ev[evid]['cnt']:
        raise Exception("no node updates for event %s" % evid)
    while True:
        if _nodes > ev[evid]['cnt']:
            break
        _inst = '%04d' % (_nodes)
        nts = UTCDateTime(ev[evid][_inst]['timest'])
        if nts >= _vts:
            break
        inst = _inst
        _nodes += 1
    print ev[evid]['cnt'], inst

    picked_stats = []
    g = pyproj.Geod(ellps='sphere')
    print "Picked stations for event %s" % evid
    for _pst in ev[evid][inst]['pickedst']:
        lon, lat, nm = _pst
        for _st in ev[evid][inst]['stations']:
            _ln, _lt, _n = _st
            az, baz, dist = g.inv(_ln, _lt, lon, lat)
            if dist < 1.0:
                nm = _n.replace('\'', '')
                break
        print "--> name: %s; lat: %.2f; lon: %.2f" % (nm, lat, lon)
        picked_stats.append(nm)
    print len(picked_stats), nnodes
    return picked_stats


def number_of_picks(reports, datdir, fout, new=True):
    if new:
        ei = {}
        p = ReportsParser()
        for _f in reports:
            p.read_reports(_f)
        correct = p.get_correct(mmin=3.5, mmax=10.0)
        for _r in correct:
            pix = get_event_info(UTCDateTime(_r[1]), os.path.basename(_r[-1]),
                                  datdir, _r[0])
            ei[_r[0]] = pix
        fh = open(fout, 'w')
        json.dump(ei, fh, indent=0)
        fh.close()
    else:
        fh = open(fout)
        ei = json.load(fh)
        fh.close()
    return ei

if __name__ == '__main__':
    pass
