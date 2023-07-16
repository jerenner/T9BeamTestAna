#!/usr/bin/python3


momentasDict = {}
runsDict = {}

# saved html page
infile = open('share/wcte-daq.html')
for xline in infile.readlines():
    line = xline[:-1]
    if not 'True' in line:
        continue
    line = line.replace('<td>','&').replace('</td>','')
    #print(line)
    tokens = line.split('&')
    #print(tokens)
    srun = tokens[1]
    smomentum = tokens[7]
    #print(srun, smomentum)
    run = int(srun)
    momentum = int(float(smomentum)*100)*10

    if run < 255:
        continue

    if abs(momentum) < 50:
        continue

    #print(run, momentum)

    try:
        n = len(momentasDict[momentum])
    except:
        momentasDict[momentum] = []
    if not run in momentasDict[momentum]:
        momentasDict[momentum].append(run)
        runsDict[run] = momentum


print('#!/usr/bin/python')
print('momentaDict = ', momentasDict)
print('runsDict = ', runsDict)


