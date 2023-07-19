#!/usr/bin/python3

from math import sqrt, pow
import ROOT
from collections import OrderedDict

momentaMeV = [200 + i*20 for i in range(1,9)]
momentaMeV.extend([500, 600, 700, 800, 900, 1000])
print(momentaMeV)

ms = OrderedDict()
ms['e'] = 0.511
ms['mu'] = 105.66
ms['pi'] = 139.57
ms['p'] = 938.3
ms['d'] = 1876.


pcols = OrderedDict()
pcols['e'] = ROOT.kRed
pcols['mu'] = ROOT.kBlue 
pcols['pi'] = ROOT.kGreen+1
pcols['p'] = ROOT.kBlack
pcols['d'] = ROOT.kGray


conv = 1.e9
#TB2022 l = 2.90 # m
#TB2023
l = 3.45 # m
c = 299792458 # m/c

print('*** The Time of Flight times for L = {} m as function of particles momenta ***'.format(l))

############################################################
def getTof(m, momentum):
    return l/c*sqrt(1.+pow(m/momentum,2))*conv

############################################################

def getTofDiff(particle1, particle2, momentum):
    m1 = ms[particle1]
    m2 = ms[particle2]
    t1 = getTof(m1, momentum)
    t2 = getTof(m2, momentum)
    return t2 - t1

############################################################

def GenerateTimes(verbose = True):

    dts = OrderedDict()
    ts = OrderedDict()
    
    for pname in ms:
        ts[pname] = []
        dts[pname] = []


    for p in momentaMeV:
        line = 'p = {} MeV: '.format(p)
        for pname in ms:
            m = ms[pname]
            t = getTof(m, p) #l/c*sqrt(1.+pow(m/p,2))*conv
            line = line + ' {}: {:2.2f}'.format(pname, t)
            ts[pname].append(1.*t)
        line = line + ' ns'
        if verbose:
            print(line)

    for pname in ts:
        if pname == 'e':
            continue
        for ip in range(0,len(ts[pname])):
            dts[pname].append(ts[pname][ip] - ts['e'][ip])

    for pname in dts:
        if pname == 'e':
            continue
        if verbose:
            print('--- dt between {} and electrons ---'.format(pname))
        for ip in range(0,len(ts[pname])):
            if verbose:
                print('p={}: dt = t_{} - t_e = {:1.2f} ns'.format(momentaMeV[ip], pname, dts[pname][ip]))

    return dts

