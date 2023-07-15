#!/usr/bin/python3

from math import sqrt, pow

from collections import OrderedDict

momentaMeV = [200 + i*20 for i in range(1,9)]
momentaMeV.extend([500, 600, 700, 800, 900, 1000])
print(momentaMeV)

ms = OrderedDict()
ms['e'] = 0.511
ms['mu'] = 135.6
ms['pi'] = 139.6
ms['p'] = 938.3
ms['d'] = 1876.

conv = 1.e9
#TB2022 l = 2.90 # m
#TB2023
l = 3.45 # m
c = 299792458 # m/c

print('*** The Time of Flight times for L = {} m as function of particles momenta ***'.format(l))

############################################################

def getTofDiff(particle1, particle2, momentum):
    m1 = ms[particle1]
    m2 = ms[particle2]
    t1 = l/c*sqrt(1.+pow(m1/momentum,2))*conv
    t2 = l/c*sqrt(1.+pow(m2/momentum,2))*conv
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
            t = l/c*sqrt(1.+pow(m/p,2))*conv
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

