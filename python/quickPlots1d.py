#!/snap/bin/pyroot

#/usr/bin/python3

# jk
# 20/09/2022
# 14.7.2023

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

from labelTools import *

cans = []
stuff = []
lines = []


def makeLine(x1, x2, y1, y2):
    line = ROOT.TLine(x1, y1, x2, y2)
    line.SetLineColor(ROOT.kGreen)
    line.SetLineWidth(2)
    line.Draw()
    return line


def PrintUsage(argv):
    print('Usage:')
    print('{} filename_plots.root [-b]'.format(argv[0]))
    print('Example:')
    print('{} output_300n_plots.root -b'.format(argv[0]))
    return

##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    pngdir = 'png_results/'
    pdfdir = 'pdf_results/'
    os.system(f'mkdir {pngdir}')
    os.system(f'mkdir {pdfdir}')

    opt2d = 'colz'
    
    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    #gBatch = True
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
            print('OK, running in batch mode')
        elif opt in ("-t", "--tag"):
            gTag = arg
            print('OK, using user-defined histograms tag for output pngs {:}'.format(gTag,) )

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    if len(argv) < 2:
        PrintUsage(argv)
        return

    ROOT.gStyle.SetOptFit(111)
    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))


    ROOT.gStyle.SetPalette(ROOT.kSolar)
    
    
    #filename = 'output_300n_plots.root'
    filename = argv[1]
    rfile = ROOT.TFile(filename, 'read')
    hbasenames = {
        'hRef_Time' : ROOT.kGreen,
        'hRef_Charge' : ROOT.kCyan,
        'hRef_Voltage' : ROOT.kMagenta,
        'hRef_nPeaks' : ROOT.kYellow,
    }
    
    nChannels = 19 # 32
    Hs = []
    Txts = []
  

    ftag = filename.split('/')[-1].replace('output_','').replace('_plots.root','')

    os.system('mkdir -p pdf png')
    
    for hbasename in hbasenames:

        hs = []
        txts = []
           
        for ich in range(0, nChannels):
            # hack just for old tofs
            #if not ( ich >= 8 and ich <= 15):
            #    continue
            hname = hbasename + str(ich)
            h = rfile.Get(hname)
            try:
                print('ok, got ', h.GetName())
            except:
                print('ERROR getting histo {}!'.format(hname))
                continue

            print('Pushing ', ich, hname)
            hs.append(h)
        Hs.append(hs)

        canname = 'WCTEJuly2023_Quick1D_{}_{}'.format(ftag, hbasename)
        canname = canname.replace('_list_root','').replace('_ntuple','')
        #can = ROOT.TCanvas(canname, canname, 0, 0, 1600, 800)
        #cans.append(can)
        #can.Divide(8,4)
        for h in hs:
            try:
                print('ok, got ', h.GetName())
            except:
                print('ERROR getting histo!')
                continue
            ich = hs.index(h)

            if ich % 8 == 0:
                idigi = ich / 8
                off = 60
                cw = 4*400 + 4*off
                ch = 2*400
                if idigi > 1:
                    cw = 2*400 + 2*off
                    ch = 400
                can = ROOT.TCanvas(canname + f'_digi{idigi}', canname + f'_digi{idigi}', 0, 0, cw, ch)
                if idigi < 2:
                    can.Divide(4,2)
                else:
                    can.Divide(3,1)
                cans.append(can)

            
            can.cd(ich % 8 + 1)
            h.SetStats(0)
            #if not 'Time' in h.GetName():
            if 'nPeaks' in h.GetName():
                ROOT.gPad.SetLogy(1)
            #h.GetYaxis().SetRangeUser(1.e-4, h.GetYaxis().GetXmax())
            h.SetFillColor(hbasenames[hbasename])
            h.SetFillStyle(1111)
            h.SetTitle(ChNames[hs.index(h)])

         
            h.Draw('hist')

      

##################################
#       plots all the canvas     #
##################################

    srun = ''
    tokens = filename.split('_')
    for token in tokens:
        if '00' in token:
            srun = token.replace('000','')
    pnote = makeMomentumLabel(srun)
    stuff.append(pnote)
    
    for can in cans:
        can.cd()
        if 'vs' in can.GetName():
            pnote.Draw()
        can.Update()
        can.Print(pngdir + can.GetName() + '.png')
        can.Print(pdfdir + can.GetName() + '.pdf')
    
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

