#!/usr/bin/env python
'''
Parse VS(EW) logfiles to extract event and pick information.
Created on May 21, 2012

@author: behry
'''
import os
from collections import defaultdict
import re
import bz2
import pyproj


def read_eewvs_log(fin, format='new', invalid=False):
    events = defaultdict(lambda: defaultdict(dict))
    picks = {}
    magcomptime = []

    # Regular expressions
    ######################

    # The following lines are written if an event has been re-evaluated taking
    # new picks and waveforms into account.
    eventpat = r'vsevent:\(quake:(?P<evid>\d+), timestamp:(?P<ts>\S+), '
    eventpat += r'nodes:(?P<nnodes>\d+), estimates:(?P<nest>\d+).*, valid:(?P<vd>\d).*'
    eventpat += r'; mag:(?P<mag>\d+.\d*).*, coord\(lat:(?P<lat>-?\d+.\d*), '
    eventpat += r'lon:(?P<lon>-?\d+.\d*)\), depth:(?P<dep>\d+\.\d*|\d+),\s+.*'
    eventpat += r'rtime:(?P<rtime>\S+), vtime:(?P<vtime>\S+),.*'
    eventpat += r'auditmagnitudeguess:(?P<mguess>\d+.\d*), '
    eventpat += r'auditmagnitudeerror:(?P<merr>\d+.\d*), '
    eventpat += r'validmagnitudeerrorratio:(?P<merrthresh>\d+.\d*), '
    eventpat += r'auditpickmissrate:(?P<pickerr>\d+.\d*|\d+), (?:azimuthal gap:(?P<azgap>\d+\.\d*), )?'
    eventpat += r'valid:(?P<vd1>\d), pgestimates:(?P<npg>\d+)'
    nodespat = r'(\d+)\s+(\S+)\s+(\S+)\s+PERL\s+(.*)'
    ppat = r'\$p->add\((-?\d+.\d*),(-?\d+.\d*),(\S*)\)'
    epat = r'\$e->add\((-?\d+.\d*),(-?\d+.\d*)\)'
    spat = r'\$s->add\((-?\d+.\d*),(-?\d+.\d*),(\S*)\)'
    wpat = r'\$w->add\((-?\d+.\d*),(-?\d+.\d*)\)'


    # This line is written if a report file is saved to disk and send by email.
    publishpat = r'publish; writing report:(\S+); \(quake:(\d+).*'

    # This line is written every time a new pick is added to an event.
    nodelist = r'update node list; quake:(?P<id>\d+), sourceid:\S+, time:(?P<ts>\S+), .*, picklist:(?P<pl>.*)'
    nodelist += r', creation time:(?P<ct>\S+)'

    # This line is written every time a pick has passed the quality criteria
    if format == 'special':
        pickpat = r'pick; linking; .*, key:(?P<pk>\d+), picktime:(?P<pt>\S+), creation time:(?P<ct>\S+), velocity.* network:(?P<net>\S+), '
    else:
        pickpat = r'pick; linking; .*, key:(?P<pk>\d+), picktime:(?P<pt>\S+), .*, acceleration:(?P<acc>\S+),.* network:(?P<net>\S+), '
    pickpat += r'station:(?P<stat>\S+),.*coord\{lat:(?P<lat>-?\d+.\d*), lon:(?P<lon>-?\d+.\d*)\}.*'

    # This line is written in regular intervals
    timeupdate = r'rtime:(\S+), vtime:(\S+), ptime:(\S+)'

    # This line is only written if the code has been run in a special modus to
    # monitor the time it takes to compute a magnitude estimate.
    comptime = r'computation time; estimate creation time: (?P<ect>\S+); '
    comptime += r'location creation time: (?P<lct>\S+); estimates:(?P<eno>\d+)'

    # This line has been added to see the reference node list that is used to
    # compute the pickmissrate
    reflist = 'Reference nodes;quake:(?P<id>\d+);(?P<nodes>.*)'

    # de-compress if necessary
    if fin.endswith('.bz2'):
        ifname = bz2.BZ2File(fin)
    else:
        ifname = open(fin)

    while True:
        line = ifname.readline()
        if not line:
            break

        match = re.search(reflist, line)
        if match:
            evid = int(match.group('id'))
            if match.group('nodes') is not None:
                if not events[evid].has_key('reflist'):
                    events[evid]['reflist'] = {}
                for _e in match.group('nodes').split(';'):
                    matchr = re.search('\$R->add\((?P<lon>-?\d+.\d*),(?P<lat>-?\d+.\d*),(?P<nm>\S*)\\)', _e)
                    if matchr:
                        lon = float(matchr.group('lon'))
                        lat = float(matchr.group('lat'))
                        statname = matchr.group('nm')
                        if statname not in events[evid]['reflist'].keys():
                            events[evid]['reflist'][statname] = (lon, lat)

        match = re.search(nodespat, line)
        if match:
            evid = int(match.group(1))
            ot = match.group(2)
            vtime = match.group(3)
            if events[evid].has_key('cnt'):
                events[evid]['cnt'] += 1
            else:
                events[evid]['cnt'] = 1
            inst = '%04d' % events[evid]['cnt']
            events[evid][inst]['origint'] = ot
            events[evid][inst]['timest'] = vtime
            events[evid][inst]['pickedst'] = []
            events[evid][inst]['stations'] = []
            events[evid][inst]['wavefront'] = []
            for _e in match.group(4).split(';'):
                matchp = re.search(ppat, _e)
                if matchp:
                    staname = matchp.group(3)
                    lon = float(matchp.group(1))
                    lat = float(matchp.group(2))
                    # print "picked station for event %d: longitude = %.3f, latitude = %.3f." % (evid, lon, lat)
                    events[evid][inst]['pickedst'].append((lon, lat, staname))
                matche = re.search(epat, _e)
                if matche:
                    lon = float(matche.group(1))
                    lat = float(matche.group(2))
                    events[evid][inst]['epic'] = (lon, lat)
                matchs = re.search(spat, _e)
                if matchs:
                    lon = float(matchs.group(1))
                    lat = float(matchs.group(2))
                    name = matchs.group(3)
                    events[evid][inst]['stations'].append((lon, lat, name))
                matchw = re.search(wpat, _e)
                if matchw:
                    lon = float(matchw.group(1))
                    lat = float(matchw.group(2))
                    events[evid][inst]['wavefront'].append((lon, lat))

        match = re.search(eventpat, line)
        if match:
            evid = int(match.group('evid'))
            ts = match.group('ts')
            nnodes = int(match.group('nnodes'))
            estno = int(match.group('nest'))
            valid = int(match.group('vd'))
            mag = float(match.group('mag'))
            lat = float(match.group('lat'))
            lon = float(match.group('lon'))
            depth = float(match.group('dep'))
            rtime = match.group('rtime')
            vtime = match.group('vtime')
            mguess = float(match.group('mguess'))
            merr = float(match.group('merr'))
            merrthresh = float(match.group('merrthresh'))
            pickerr = float(match.group('pickerr'))
            valid1 = int(match.group('vd1'))
            npg = int(match.group('npg'))
            if match.group('azgap') is not None:
                azgap = float(match.group('azgap'))
            else:
                azgap = -1
            if valid1 != 1 and not invalid:
                pass
            else:
                if not events.has_key(evid) or \
                not events[evid].has_key('estimates'):
                    events[evid]['estimates'] = []
                events[evid]['estimates'].append((estno, nnodes, mag, lat, lon, depth, rtime, vtime, ts, mguess, merr, merrthresh, pickerr, azgap))

        match = re.search(publishpat, line)
        if match:
            repfile = match.group(1)
            evid = int(match.group(2))
            events[evid]['report'] = repfile

        match = re.search(pickpat, line)
        if match:
            pickkey = int(match.group('pk'))
            ptime = match.group('pt').split('(')[0]
            if format == 'special':
                ctime = match.group('ct')
            else:
                ctime = 0.0
            net = match.group('net').split('(')[0]
            stat = match.group('stat').split('(')[0]
            lat = float(match.group('lat'))
            lon = float(match.group('lon'))
            acc = float(match.group('acc'))
            if pickkey not in picks.keys():
                picks[pickkey] = []
            # The pick ID number starts at 1000 and is reset
            # if eewvs is restarted. Therefore in one log file
            # several picks with the same ID can occur (see EewVsPick.cc).
            picks[pickkey].append((ptime, ctime, net, stat, lat, lon, acc))

        match = re.search(nodelist, line)
        if match:
            eventkey = int(match.group('id'))
            ot = match.group('ts').split('(')[0]
            try:
                picklist = map(int, match.group('pl').split())
            except Exception, e:
                print e
                print line
                print match.group('pl')
            if not events[eventkey].has_key('nodelist'):
                events[eventkey]['nodelist'] = []
            if format == 'new':
                try:
                    ct = match.group(4)
                except:
                    print line
            else:
                ct = 0.0
            events[eventkey]['nodelist'].append((ot, ct, picklist))

        match = re.search(timeupdate, line)
        if match:
            recenttime = match.group(1).split('@')[0]
            vstime = match.group(2).split('@')[0]
            picktime = match.group(3).split('@')[0]

        match = re.search(comptime, line)
        if match:
            try:
                ect = match.group('ect')
                lct = match.group('lct')
                eno = int(match.group('eno'))
                if eno == 1:
                    print eno, ect, lct, ect - lct
                    magcomptime.append(ect - lct)
            except:
                print line

    ifname.close()
    return events, picks, magcomptime


def find_delayed_stations(fin):
    """
    Find all entries reporting unavailable waveforms.
    """
    pattern = r'datagram not available;.*, picktime:(\S+),.*network:(\S+), station:(\S+),.*'
    delays = defaultdict(list)
    if fin.endswith('.bz2'):
        ifname = bz2.BZ2File(fin)
    else:
        ifname = open(fin)
    for line in ifname.readlines():
        match = re.search(pattern, line)
        if match:
            ptime = match.group(1).split('(')[0]
            net = match.group(2).split('(')[0]
            stat = match.group(3).split('(')[0]
            delays[stat].append(ptime)
    ifname.close()
    return delays


def summary(events, picks, eid, allstats=False, updates=False, estimates=False,
            nodes=False):
    """
    Print some of the information gathered from the logfiles.
    """
    print "Summary for event {:d}".format(eid)
    g = pyproj.Geod(ellps='WGS84')
    if updates:
        for i in xrange(events[eid]['cnt']):
            inst = '%04d' % (i + 1)
            ot = events[eid][inst]['origint']
            vtime = events[eid][inst]['timest']
            ps = events[eid][inst]['pickedst']
            stats = events[eid][inst]['stations']
            print "Update: {:s}".format(inst)
            print "Creation time: {:s}; origin time: {:s}".format(vtime, ot)
            print "Picked stations:"
            for st in ps:
                lon, lat, staname = st
                for _st in stats:
                    _ln, _lt, _n = _st
                    az, baz, dist = g.inv(_ln, _lt, lon, lat)
                    if dist < 1.0:
                        nm = _n.replace('\'', '')
                        break
                print "-->{:s} ({:.2f}, {:.2f})".format(nm, lon, lat)

    if nodes:
        for _nds in events[eid]['nodelist']:
            ot, ct, picklist = _nds
            print "Pick ids: ", " ".join(map(str, picklist))

    if allstats:
        print "All stations: "
        line = ''
        for _st in stats:
            _ln, _lt, _n = _st
            line += '{:s}; '.format(_n)
        print line

    if estimates:
        hdrentries = ('#', 'mag', 'lat', 'lon', 'dep', 'creation time',
                      'virtual time', 'origin time', 'mag est', 'mag diff',
                      'err thresh', 'pck err', 'nst')
        hdr = "{:^3s}{:^5s}{:^8s}{:^8s}{:^8s}{:^26s}{:^26s}{:^26s}{:^9s}{:^10s}"
        hdr += "{:^12s}{:^9s}{:^5s}"
        print hdr.format(*hdrentries)
        print "{0:-^155}".format("")
        for _est in events[eid]['estimates']:
            l = "{0:3d}{2:5.2f}{3:8.2f}{4:8.2f}{5:8.2f}"
            l += "{6:>26s}{7:>26s}{8:>26s}{9:9.2f}{10:10.3f}"
            l += "{11:12.3f}{12:9.3f}{1:4d}"
            print l.format(*_est)


def showpicks(picks):
    """
    Show the picks that passed the quality criteria.
    """
    hdrentries = ('ID', 'pick time', 'net', 'sta', 'lat', 'lon', 'acc (cm/s/s)')
    hdr = "{:^8s}{:^25s}{:^5s}{:^5s}{:^8s}{:^8s}{:^12s}"
    print hdr.format(*hdrentries)
    print "{0:-^71}".format("")
    for _k in picks.keys():
        for pk in picks[_k]:
            pt, ct, net, stat, lat, lon, acc = pk
            l = "{:8d}{:>25s}{:>5s}{:>5s}{:8.2f}{:8.2f}{:12.4e}"
            print l.format(_k, pt, net, stat, lat, lon, acc)


if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", help="Name of the report file to look for.")
    parser.add_argument("--evid", help="Event ID of the event you are interested in.")
    parser.add_argument("--invalid", help="Include invalid estimates in addition to valid estimates.",
                        action='store_true')
    parser.add_argument("--updates", help="Show station updates.",
                        action='store_true')
    parser.add_argument("--estimates", help="Show the magnitude and location estimates.",
                        action='store_true')
    parser.add_argument("--pickids", help="Show the IDs of the picks that contributed to an event.",
                        action='store_true')
    parser.add_argument("--picks", help="Show the picks received by VS(EW).",
                        action='store_true')
    parser.add_argument('fn', metavar='Log-file', type=str,
                        help='The path to the log-file to analyse')
    args = parser.parse_args()
    events, picks, magcomptime = read_eewvs_log(args.fn, invalid=args.invalid,
                                                format='new')
    if args.picks:
        showpicks(picks)
        sys.exit(0)
    if args.report:
        for _k in sorted(events.keys()):
            try:
                if events[_k]['report'] == args.report:
                    summary(events, picks, _k, updates=args.updates,
                            estimates=args.estimates, nodes=args.pickids)
            except Exception, e:
                print e
                continue
    elif args.evid:
        eid = int(args.evid)
        try:
            summary(events, picks, eid, updates=args.updates,
                    estimates=args.estimates, nodes=args.pickids)
        except Exception, e:
            print e
    else:
        for _k in sorted(events.keys()):
            try:
                summary(events, picks, _k, updates=args.updates,
                        estimates=args.estimates, nodes=args.pickids)
            except Exception, e:
                print e
                continue

