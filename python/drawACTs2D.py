#!/usr/bin/python3
# 20/09/2022

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

def PrintUsage(argv):
    print('Usage:')
    print('{} filename_plots.root A/C'.format(argv[0]))
    print('Example:')
    print('{} output_300n_plots.root A -b'.format(argv[0]))
    return

##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    #gBatch = True
    gBatch = False
    gTag=''
    print(argv[1:])
    try:
        # options that require an argument should be followed by a colon (:).
        opts, args = getopt.getopt(argv[3:], 'hbt:', ['help','batch','tag='])
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
            print('OK, running in batch mode')
        elif opt in ("-t", "--tag"):
            gTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gTag,) )

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    if len(argv) < 2:
        PrintUsage(argv)
        return

    basetag = 'charges'
    
    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetPadLeftMargin(0.20)
    
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    #ROOT.gStyle.SetPalette(ROOT.kSolar)
    
    #filename = 'output_300n_plots.root'
    filename = argv[1]
    rfile = ROOT.TFile(filename, 'read')
    hnames = [ 'hRef_ACT2CACT1C', 'hRef_ACT3CACT2C', 'hRef_ACT1CACT3C']
    hs = []
    txts = []


    ftag = filename.replace('output_','').replace('_plots.root','')
    canname = 'WCTEJuly2022_ACT2D_{}_{}'.format(basetag, ftag)
    can = ROOT.TCanvas(canname, canname, 0, 0, 1500, 450)
    can.Divide(3,1)
    cans.append(can)

    
    for hname in hnames:
        h = rfile.Get(hname)
        hs.append(h)
    for h in hs:
        hnameTag = h.GetName().replace('hRef_','')
        can.cd(hs.index(h)+1)
        h.SetStats(0)
        #ROOT.gPad.SetLogy(1)
        #h.GetYaxis().SetRangeUser(1.e-4, h.GetYaxis().GetXmax())
        h.Draw('colz')
        ROOT.gPad.SetLogz(1)
        txt= '{} {} p={} MeV/c'.format(hnameTag, basetag, ftag.replace('n','').replace('p',''))
        if 'p' in ftag:
            txt = txt + ' (Pos)'
        else:
            txt = txt + ' (Neg)'
        text = ROOT.TLatex(0.120, 0.93, txt)
        text.SetTextSize(0.04)
        text.SetNDC()
        text.Draw()
        txts.append(text)
       

    can.Update()
    can.Print(can.GetName() + '.png')
    can.Print(can.GetName() + '.pdf')
    if not gBatch:
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

