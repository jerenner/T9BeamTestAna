#!/snap/bin/pyroot
# was: #!/usr/bin/python3
# Pá 14. července 2023, 18:39:38 CEST

#from __future__ import print_function

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

def PrintUsage(argv):
    print('Usage:')
    print('{} histos/output_list_root_run_XYZ_plots.root [-b]'.format(argv[0]))
    return


def Fit(h, tag, ct, w, t1, t2, peaksf = 1.):
    fname = 'fit{}'.format(tag)
    hname = h.GetName()
    fit = ROOT.TF1(fname, '[0]*exp(-(x-[1])^2/(2*[2]^2))', t1, t2)
    #ampl = h.GetMaximum() / peaksf
    ampl = h.GetBinContent(h.FindBin(ct))
    print('Amplitude initially ', ampl)
    fit.SetParameters(ampl, ct, w)
    h.Fit(fname, 'q', '')
    mean = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    #print(f'1) {hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
    sf = 2.
    h.Fit(fname, 'q', '', mean - sf*sigma, mean + sf*sigma)
    mean = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    print(f'{hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
    fit.SetLineColor(h.GetLineColor())
    fit.SetLineStyle(2)
    fit.Draw('same')
    return fit






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

    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))


    if len(argv) < 2:
        PrintUsage(argv)
        return


    fileName = argv[1]
    inFile = ROOT.TFile(fileName, "READ")
        
    hTOFAll = inFile.Get("hTOFAll")
    hTOFEl = inFile.Get("hTOFEl")


    canname = 'FitToF'
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    cans.append(can)


    #can.Divide(2,2)
    
    ROOT.gStyle.SetOptStat(0)
    
    hTOFAll.SetLineColor(ROOT.kBlack)
    hTOFAll.SetLineWidth(2)
    hTOFAll.Draw('hist')

    # def Fit(h, tag, ct, w, t1, t2, peaksf = 1.):

    off = 0.
    if 'uncalibrated' in inFile.GetName():
        off = 3.
    fite = Fit(hTOFAll, '_el', 11.7 + off, 1., 11. + off, 13. + off, 1.)
    fitp = Fit(hTOFAll, '_p',  16. + off, 1., 14 + off,  17.5 + off, 1.)
    #fitd = Fit(hTOFAll, '_d',  25., 1., 24, 26, 1.)

    tdif_e_p = fitp.GetParameter(1) - fite.GetParameter(1)
    tex = ROOT.TLatex(0.7, 0.8, 't_{p}-t_{e}=' + '{:1.2f}'.format(tdif_e_p))
    tex.SetNDC()
    tex.Draw()
    stuff.append(tex)

    #tdif_e_d = fitd.GetParameter(1) - fite.GetParameter(1)
    #tex = ROOT.TLatex(0.7, 0.7, 't_{d}-t_{e}=' + '{:1.2f}'.format(tdif_e_d))
    #tex.SetNDC()
    #tex.Draw()
    #stuff.append(tex)


    
    hTOFEl.SetLineColor(ROOT.kRed)
    hTOFEl.SetLineWidth(2)
    hTOFEl.Draw('hist same')

    ROOT.gPad.SetLogy(1)
    ROOT.gPad.Update()
    
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

