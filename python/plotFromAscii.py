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
    
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))
    
    xafilenames = os.popen('ls ascii*.txt {}'.format(addgrep)).readlines()


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
        pstr = alines[0].split()[1] + alines[0].split()[0]
        
        for line in alines:
            for part in parts:
                if part in line.split()[0]:
                    Ns[part].append([pstr, float(line.split()[1]) ] )
    print(Ns)

    grs = {}
    for part in parts:
        gr = ROOT.TGraph()
        gr.SetName('gr_' + part)
        grs[part] = gr

    for part in parts:
        gr = grs[part]
        for pair in Ns[part]:
            j = gr.GetN()
            pstr = pair[0].replace('-','')
            ipstr = pstr.replace('p','').replace('n','')
            momentum = abs(int(ipstr))
            n = pair[1] / data_runs.Runs[pstr][1]
            gr.SetPoint(j, momentum, n)

    canname = 'MultiPerSpill_TBJuly2022_' + signTag
    can = ROOT.TCanvas(canname, canname, 100, 100, 1100, 800)
    cans.append(can)

    pmin = 150.
    pmax = 400.
    h2 = ROOT.TH2D("tmp", "tmp;p [MeV/c];N / spill", 100, pmin, pmax, 100, 0, 50.e1)
    h2.SetStats(0)
    h2.Draw()

    leg = ROOT.TLegend(0.55, 0.7, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetHeader('WCTE TB July 2022')
    opt = 'PL'
    for part in grs:
        gr = grs[part]
        gr.SetMarkerColor(parts[part])
        gr.SetLineColor(parts[part])
        gr.SetMarkerStyle(20)
        gr.Draw(opt)
        leg.AddEntry(gr, 'N_{' + part + '} / spill (' + signTag + ')', 'PL')
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

