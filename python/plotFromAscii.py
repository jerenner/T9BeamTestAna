#!/usr/bin/python3
from collections import OrderedDict
import data_runs
import pytools

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

import ctypes

cans = []
stuff = []


def printGr(gr):
    x = ctypes.c_double(0.)
    y = ctypes.c_double(0.)
    for ip in range(0, gr.GetN()):
        gr.GetPoint(ip, x,y)
        print(ip, x, y)

##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):

    # CONVERSION CONSTANTS
    spillSep = 40. # s
    secsInDay = 24*3600
    tsf = secsInDay / spillSep

    ROOT.gStyle.SetPadLeftMargin(0.15)
    
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

    if len(argv) < 3:
        print('Usage: need to specify n or p for negative or positive momenta to plot the mu and pi yields over momenta!')
        print('{} n/p low/high/p'.format(argv[0]))
        exit(1)
    

    addgrep = ''
    signTag = 'all'
    print(argv)
    if argv[1] == 'n' or argv[1] == '_n' or argv[1] == 'neg' or argv[1] == 'Neg' or argv[1] == '_neg' or argv[1] == '_Neg':
        addgrep = ' | grep n_'
        signTag = 'Neg'
    if argv[1] == 'p' or argv[1] == '_p' or argv[1] == 'pos' or argv[1] == 'Pos' or argv[1] == '_pos' or argv[1] == '_Pos':
        addgrep = ' | grep p_'
        signTag = 'Pos'

    isProtonFit = False
    Runs = pytools.getRuns(argv[2])
    twoComponentFit = argv[2] == 'high'
    
    if argv[2] == 'p' and signTag != 'p':
        print('Asked for proton runs, so switching to assuming positive momenta!')
        signTag = 'Pos'
        isProtonFit = True

    print('Configured as: addreg="{}", sigTag={}'.format(addgrep, signTag))
             
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    egrep = ' | egrep "'
    irun = 0
    for run in Runs:
        egrep = egrep + run
        if irun < len(Runs) - 1:
            egrep = egrep + '|'
        irun = irun + 1
    egrep = egrep + '"'
    print(egrep)
    xafilenames = os.popen('ls ascii*.txt {} {}'.format(egrep,addgrep)).readlines()
    print(xafilenames)
    if '1000' in xafilenames[0]:
        xafilenamestmp = xafilenames[1:]
        xafilenamestmp.append(xafilenames[0])
        xafilenames = xafilenamestmp
    print(xafilenames)
    parts = {'e' : ROOT.kRed,
             'mu' : ROOT.kBlue,
             'pi' : ROOT.kGreen+2}
    if twoComponentFit:
        parts = {'e' : ROOT.kRed,
                 'mupi' : ROOT.kMagenta
        }
    
    if isProtonFit:
        parts =  {'e' : ROOT.kRed,
                  'proton' : ROOT.kBlack}
    
    Ns = OrderedDict()
    Effs = OrderedDict()
    for part in parts:
        Ns[part] = []
        Effs[part] = []
    for xafilename in xafilenames:
        afilename = xafilename[:-1]
        print(afilename)
        afile = open(afilename, 'r')
        alines = []
        
        for aline in afile.readlines():
            alines.append(aline)
        signStr = 'p'
        if int(alines[0].split()[1]) < 0:
            signStr = 'n'
        pstr = alines[0].split()[1] + signStr
        #print(pstr)

        iline = -1
        for line in alines:
            print(line)
            iline = iline + 1
            print('Processing particle {}'.format(part))
            for part in parts:
                if line.split()[0] == 'N_' + part:
                    Ns[part].append([pstr, float(line.split()[1]), float(line.split()[2])] )
                if line.split()[0] == 'eff_' + part:
                    Effs[part].append([pstr, float(line.split()[1]), float(line.split()[2])] )
    print(Ns)
    print(Effs)

    # Expected particle counts:

    grs = OrderedDict()
    for part in parts:
        gr = ROOT.TGraphErrors()
        gr.SetName('gr_' + part)
        grs[part] = gr

    for part in parts:
        gr = grs[part]
        for data in Ns[part]:
            j = gr.GetN()
            pstr = data[0].replace('-','')
            ipstr = pstr.replace('p','').replace('n','')
            momentum = abs(int(ipstr))
            nspills  = Runs[pstr][1]
            print('p={} spills={}'.format(pstr, nspills))
            n = data[1] / nspills * tsf
            err = data[2] / nspills * tsf
            print(part, n)
            gr.SetPoint(j, momentum, n)
            gr.SetPointError(j, 0, err)

    print(grs)
    canname = 'MultiPerSpill_TBJuly2022_' + signTag + "_" + argv[2]
    can = ROOT.TCanvas(canname, canname, 100, 100, 1100, 800)
    cans.append(can)

    pmin = 180.
    pmax = 300.
    y0 = 0
    y1 = 20 * tsf

    if twoComponentFit:
        pmin = 290.
        pmax = 370.
        y1 = 80 * tsf
    
    if isProtonFit:
        pmin = 300.
        pmax = 1100.
        y1 = 1.5*300 * tsf

    h2 = ROOT.TH2D("tmp", "tmp;|p| [MeV/c];N / day", 100, pmin, pmax, 100, y0, y1)
    h2.SetStats(0)
    h2.Draw()

    ROOT.gPad.SetGridy(1)
    
    leg = ROOT.TLegend(0.55, 0.7, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetHeader('WCTE TB July 2022')
    opt = 'PLe1'
    for part in grs:
        gr = grs[part]
        gr.SetMarkerColor(parts[part])
        gr.SetLineColor(parts[part])
        gr.SetLineStyle(1)
        gr.SetLineWidth(2)
        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(1.2)
        if part != 'e':
            print('Drawing graph for {}'.format(part))
            gr.Draw(opt)
            #printGr(gr)
            texX = part
            if part == 'pi' or part == 'mu':
                texX = '#' + part
            if part == 'mupi':
                texX = '#mu+#pi'
            leg.AddEntry(gr, 'N_{' + texX + '} / day (' + signTag + ')', 'PL')
    leg.Draw()

    
    can.Update()
    can.Print(can.GetName() + '.png')
    can.Print(can.GetName() + '.pdf')


    # Effs:

    grsEff = OrderedDict()
    for part in parts:
        if part == 'e':
            continue
        gr = ROOT.TGraphErrors()
        gr.SetName('grEff_' + part)
        grsEff[part] = gr

    for part in parts:
        if part == 'e':
            continue
        gr = grsEff[part]

        for data in Effs[part]:
            #print('data ', data)
            j = gr.GetN()
            pstr = data[0].replace('-','')
            ipstr = pstr.replace('p','').replace('n','')
            momentum = abs(int(ipstr))
            eff = data[1]
            err = data[2]
            print(part, eff)
            gr.SetPoint(j, momentum, eff)
            gr.SetPointError(j, 0, err)

    print(grsEff)
    canname = 'Effs_TBJuly2022_' + signTag + "_" + argv[2]
    canEff = ROOT.TCanvas(canname, canname, 100, 100, 1100, 800)
    cans.append(canEff)

    pmin = 180.
    pmax = 300.
    y0 = 0
    y1 = 1.4

    if twoComponentFit:
        pmin = 290.
        pmax = 370.
    
    if isProtonFit:
        pmin = 300.
        pmax = 1100.

    h2Eff = ROOT.TH2D("tmpEff", "tmpEff;|p| [MeV/c];Efficiency", 100, pmin, pmax, 100, y0, y1)
    h2Eff.SetStats(0)
    h2Eff.Draw()

    ROOT.gPad.SetGridy(1)
    
    legEff = ROOT.TLegend(0.55, 0.7, 0.88, 0.88)
    legEff.SetBorderSize(0)
    legEff.SetHeader('WCTE TB July 2022')
    opt = 'PLe1'
    for part in grsEff:
        print('Drawing eff for {}'.format(part))
        gr = grsEff[part]
        gr.SetMarkerColor(parts[part])
        gr.SetLineColor(parts[part])
        gr.SetLineStyle(1)
        gr.SetLineWidth(2)
        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(1.2)
        if part != 'e':
            print('Drawing graph for {}'.format(part))
            gr.Draw(opt)
            #printGr(gr)
            texX = part
            if part == 'pi' or part == 'mu':
                texX = '#' + part
            if part == 'mupi':
                texX = '#mu+#pi'
            legEff.AddEntry(gr, '#epsilon_{' + texX + '}^{' + signTag + '}', 'PL')
    legEff.Draw()

    
    canEff.Update()
    canEff.Print(canEff.GetName() + '.png')
    canEff.Print(canEff.GetName() + '.pdf')
    
    
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

