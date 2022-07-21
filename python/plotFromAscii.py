#!/usr/bin/python3

import data_runs

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []


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

    if len(argv) < 2:
        print('Usage: need to specify n or p for negative or positive momenta to plot the mu and pi yields over momenta!')
        print('{} n/p'.format(argv[0]))
        exit(1)
    

    addgrep = ''
    signTag = 'all'
    print(argv)
    if len(argv) > 1:
         if argv[1] == 'n' or argv[1] == '_n' or argv[1] == 'neg' or argv[1] == 'Neg' or argv[1] == '_neg' or argv[1] == '_Neg':
             addgrep = ' | grep n_'
             signTag = 'Neg'
         if argv[1] == 'p' or argv[1] == '_p' or argv[1] == 'pos' or argv[1] == 'Pos' or argv[1] == '_pos' or argv[1] == '_Pos':
             addgrep = ' | grep p_'
             signTag = 'Pos'
             
    print('Configured as: addreg="{}", sigTag={}'.format(addgrep, signTag))
             
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))
    
    xafilenames = os.popen('ls ascii*.txt {}'.format(addgrep)).readlines()
    print(xafilenames)

    parts = {'e' : ROOT.kRed,
             'mu' : ROOT.kBlue,
             'pi' : ROOT.kGreen+2}
    Ns = {}
    for part in parts:
        Ns[part] = []
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
        
        for line in alines:
            for part in parts:
                if part in line.split()[0]:
                    Ns[part].append([pstr, float(line.split()[1]), float(line.split()[2])] )
    print(Ns)

    grs = {}
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
            nspills  = data_runs.Runs[pstr][1]
            print('p={} spills={}'.format(pstr, nspills))
            n = data[1] / nspills * tsf
            err = data[2] / nspills * tsf
            gr.SetPoint(j, momentum, n)
            gr.SetPointError(j, 0, err)

    canname = 'MultiPerSpill_TBJuly2022_' + signTag
    can = ROOT.TCanvas(canname, canname, 100, 100, 1100, 800)
    cans.append(can)

    pmin = 180.
    pmax = 300.
    y0 = 0
  
    y1 = 20 * tsf
    h2 = ROOT.TH2D("tmp", "tmp;p [MeV/c];N / day", 100, pmin, pmax, 100, y0, y1)
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
            gr.Draw(opt)
            texX = ''
            if part == 'pi' or part == 'mu':
                texX = '#'
            leg.AddEntry(gr, 'N_{' + texX + part + '} / spill (' + signTag + ')', 'PL')
    leg.Draw()

    
    can.Update()
    can.Print(can.GetName() + '.png')
    can.Print(can.GetName() + '.pdf')
    
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

