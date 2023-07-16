#!/usr/bin/python3

from math import sqrt, pow

from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["science", "notebook", "grid"]) #for pretty plot

momentaMeV = [100 + i*10 for i in range(1,100)]
#momentaMeV.extend([500, 600, 700, 800, 900, 940, 1000])
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


def GenerateTimes():

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
        print(line)

    for pname in ts:
        if pname == 'e':
            continue
        for ip in range(0,len(ts[pname])):
            dts[pname].append(ts[pname][ip] - ts['e'][ip])

    for pname in dts:
        if pname == 'e':
            continue
        print('--- dt between {} and electrons ---'.format(pname))
        for ip in range(0,len(ts[pname])):
            print('p={}: dt = t_{} - t_e = {:1.2f} ns'.format(momentaMeV[ip], pname, dts[pname][ip]))

    return dts, np.array(ts['e'])

def TofSeparationToBeam(tofResolution = 0):

    p_array = []
    for i in range(len(tofDiff)):
        #11.51 is the electron hit time
        m =  ms['p']
        p = m/sqrt(pow((tofDiff[i]+tofResolution+t_electron[i]) * c / (conv * l), 2) - 1)
        p_array.append(p)
    #= m/sqrt(((tofDiff+tofResolution+t_electron) * c / (conv * l)) ** 2 - 1)

        #print(p)

    return np.array(p_array)






diff_to_e, t_electron = GenerateTimes()
print(diff_to_e['p'])

tofDiff = np.array(diff_to_e['p'])
predMomentaMeV = TofSeparationToBeam()
#print(np.array(predMomentaMeV) - np.array(momentaMeV))

momentaMeV = np.array(momentaMeV)

plt.figure()
for offset in np.arange(-0.75, 1, 0.25):
    plt.plot(momentaMeV, (TofSeparationToBeam(offset)-momentaMeV)/momentaMeV, label = 'TOF difference offset = %.2fns'%offset)
plt.legend()
plt.xlabel("True beam momentum")
plt.ylabel("Estimated - True beam / True beam momentum")
plt.title("Impact of TOF uncertainty on beam momentum estimation - using e-p separation")

plt.figure()
plt.plot(momentaMeV, diff_to_e['p'], label = 'dt e-p')
plt.plot(momentaMeV, diff_to_e['mu'], label = 'dt e-mu')
plt.plot(momentaMeV, np.array(diff_to_e['pi'])-np.array(diff_to_e['mu']), label = 'dt pi-mu')
#plt.semilogy()
plt.xlabel("Beam momentum (MeV/c)")
plt.ylabel("TOF difference (ns)")
plt.legend()
# plt.grid()
plt.title("(theoretical) Difference between TOF of particles @ beam momentum 100-1000 MeV/c")

plt.show()
