#!/snap/bin/pyroot
# was: #!/usr/bin/python3
# JK St 26. července 2023, 16:45:25 CEST

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from collections import OrderedDict

MeV = 1.

lst = {-1 : 2, 1: 1}
mst = {-1 : 21, 1: 8}
cols = { 'e' : ROOT.kRed, 'mu' : ROOT.kBlue, 'pi' : ROOT.kGreen+2, 'p' : ROOT.kBlack, 'mu+pi' : ROOT.kMagenta}
psymb = { 'e' : 'e', 'mu' : '#mu', 'pi' : '#pi', 'p' : 'p', 'mu+pi' : '#mu+#pi'}


####################################################################################

def getFittedResults(infilename):
    infile = open(infilename)
    results = {}
    for xline in infile.readlines():
        line = xline[:-1]
        if not 'Momentum' in line:
            continue
        tokens = line.split(',')
        sp = tokens[0].replace('Momentum','').replace('MeV/c','').replace(' ','')
        p = int(sp)

        target = tokens[1].split(':')[1].replace(' ','')
        itarget = -1
        if target == 'Alu':
            target = 'T3 Al'
            itarget = 3
        elif target == 'Tun':
            target = 'T1 Be+W'
            itarget = 1

        nSpills = float(tokens[6].split(':')[1].replace(' ',''))
        print(f'nSpills: {nSpills}')
        nmus = tokens[4].split(':')
        nmu = int(nmus[1]) / nSpills
        enmu = int(nmus[2]) / nSpills
        Nmu = [nmu, enmu]

        npis = tokens[5].split(':')
        npi = int(npis[1]) / nSpills
        enpi = int(npis[2]) / nSpills
        Npi = [npi, enpi]

        result = {}
        result['mu'] = Nmu
        result['pi'] = Npi

        print(result)
        results[target] = result

    return p,results
    

cans = []
stuff = []

####################################################################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    gBatch = False
    gTag=''
    print(argv[1:])
    try:
        # options that require an argument should be followed by a colon (:).
        opts, args = getopt.getopt(argv[2:], 'hbt:', ['help','batch','tag='])

        print('Got options:')
        print(opts)
        print(args)
    except getopt.GetoptError:
        print('Parsing...')
        print ('Command line argument error!')
        print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]]'.format(argv[0]))
        sys.exit(2)
    for opt,arg in opts:
        print('Processing command line option {} {}'.format(opt,arg))
        if opt == '-h':
            print('{:} [ -h -b --batch -tTag --tag="MyCoolTag"]'.format(argv[0]))
            sys.exit()
        elif opt in ("-b", "--batch"):
            gBatch = True
        elif opt in ("-t", "--tag"):
            gTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gTag,) )

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    ROOT.gStyle.SetOptTitle(0)
        
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    results = OrderedDict()

    indir= 'fitres/'
    xinfilenames = os.popen(f'cd {indir} ; ls fitres*.txt').readlines()
    infilenames = []
    for xfile in xinfilenames:
        infilenames.append(indir + xfile[:-1])
    # HACK
    #infilenames = ['fitres/fitres_220.0.txt']
    print(infilenames)
    
    Results = {}
    for infilename in infilenames:
        p,result = getFittedResults(infilename)
        results[p] = result

        print(p, result)
        Results[p] = result

    Grs = {}
    for p,results in Results.items():
        for target, ndict in results.items():
            try:
                ngr = len(Grs[target])
            except:
                Grs[target] = {}
                
            for part, nums in ndict.items():
                try:
                    ngr = len(Grs[target][part])
                except:
                    Grs[target][part] = {}    
                sgnp = p / abs(p)
                try:
                    ngr = Grs[target][part][sgnp].GetN()
                except:
                    Grs[target][part][sgnp] = ROOT.TGraphErrors()
                    print(f'Making new TGraph for {target} and {part}, of color {cols[part]} and momentum {p} signum {sgnp}')
                    Grs[target][part][sgnp].SetLineColor(cols[part])
                    Grs[target][part][sgnp].SetMarkerColor(cols[part])
                    Grs[target][part][sgnp].SetMarkerStyle(mst[sgnp])
                    Grs[target][part][sgnp].SetMarkerSize(1)
                    Grs[target][part][sgnp].SetLineStyle(lst[sgnp]) # by signum
                    Grs[target][part][sgnp].SetLineWidth(2)
                ip = Grs[target][part][sgnp].GetN()
                Grs[target][part][sgnp].SetPoint(ip, abs(p), nums[0])
                Grs[target][part][sgnp].SetPointError(ip, 0, nums[1])
        
    print(Grs)
    canname = 'can_partsPerSpill'
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 600)
    can.Divide(2,1)
    cans.append(can)


    pmin = 180.*MeV
    pmax = 400.*MeV
    y0, y1 = 0., 15.
    
    h2 = ROOT.TH2D("tmp", "tmp;|p| [MeV/c];N / spill", 100, pmin, pmax, 100, y0, y1)
    h2.SetStats(0)
    stuff.append(h2)

    legs = {}
    txts = {}
    
    ican = 1
    for target,PartsGr in Grs.items():
        txt = ROOT.TLatex(0.12, 0.92, f'Target: {target}')
        txt.SetNDC()
        legs[target] = ROOT.TLegend(0.16, 0.75, 0.75, 0.88)
        leg = legs[target]
        leg.SetNColumns(2)
        leg.SetBorderSize(0)
        print(f'...plotting target {target}')
        can.cd(ican)
        ROOT.gPad.SetGridy(1)
        ROOT.gPad.SetGridx(1)
        h2.DrawCopy()
        for part, partsGr in PartsGr.items():
            for sgnp,gr in partsGr.items():
                sgnstr = '+'
                if sgnp < 0:
                    sgnstr = '-'
                print(f'  ...plotting {part}')
                gr.Draw('PL')
                leg.AddEntry(gr, f'{psymb[part]} {sgnstr}', 'PL')
        leg.Draw()
        txt.Draw()
        stuff.append([txt, leg])
        ican = ican + 1

    


    can.Update()
    ROOT.gApplication.Run()
    return

###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################

