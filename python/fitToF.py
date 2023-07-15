#!/snap/bin/pyroot
# was: #!/usr/bin/python3
# Pá 14. července 2023, 18:39:38 CEST

#from __future__ import print_function

from data_runs import *
from tofUtil import *

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt

cans = []
stuff = []

##########################################

def PrintUsage(argv):
    print('Usage:')
    print('{} histos/output_list_root_run_XYZ_plots.root [-b]'.format(argv[0]))
    return

##########################################

def Fit(h, tag, ct, w, t1, t2, peaksf = 1.):
    fname = 'fit{}'.format(tag)
    hname = h.GetName()
    fit = ROOT.TF1(fname, '[0]*exp(-(x-[1])^2/(2*[2]^2))', t1, t2)
    fit.SetLineColor(ROOT.kBlack)
    fit.SetLineStyle(2)
    #ampl = h.GetMaximum() / peaksf
    ampl = h.GetBinContent(h.FindBin(ct)) / peaksf
    print('Amplitude initially ', ampl)
    fit.SetParameters(ampl, ct, w)
    fit.SetParLimits(2, 0., 0.8)
    for ip in range(0, fit.GetNpar()):
        print(fit.GetParameter(ip))
    #prefit = fit.DrawCopy('same')
    #stuff.append(prefit)
    h.Fit(fname, '', '')
    mean = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    #print(f'1) {hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
    sf = 2.
    h.Fit(fname, '', '', mean - sf*sigma, mean + sf*sigma)
    mean = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    print(f'{hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
    fit.SetLineColor(h.GetLineColor())
    fit.SetLineStyle(2)
    fit.Draw('same')
    return fit


##########################################



##########################################
##########################################
##########################################

# https://www.tutorialspoint.com/python/python_command_line_arguments.htm

def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    pngdir = 'png_results/'
    pdfdir = 'pdf_results/'
    os.system(f'mkdir {pngdir}')
    os.system(f'mkdir {pdfdir}')
    
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
    srun = fileName.split('/')[-1].replace('output_list_root_run_000','').replace('_plots.root','')
    momentum = getMomentum(srun)
    print(f'Assuming run {srun} and momentum {momentum}')

    suff = ''
    if abs(momentum) < 500:
        suff = 'Low'
        print('Low momentum run, looking at zoomed version of tof histos!')
    hTOFAll = inFile.Get("hTOFAll" + suff)
    hTOFEl = inFile.Get("hTOFEl" + suff)


    canname = 'FitToF'
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    cans.append(can)


    #can.Divide(2,2)
    
    ROOT.gStyle.SetOptStat(0)
    
    hTOFAll.SetLineColor(ROOT.kBlack)
    hTOFAll.SetMarkerColor(hTOFAll.GetLineColor())
    hTOFAll.SetMarkerStyle(20)
    hTOFAll.SetMarkerSize(1)
    hTOFAll.SetLineWidth(1)
    hTOFAll.Draw('e1')

    # def Fit(h, tag, ct, w, t1, t2, peaksf = 1.):

    off = 0.
    if 'uncalibrated' in inFile.GetName():
        off = 3.
    width = 0.2
    fite = Fit(hTOFAll, '_el', 11.7 + off, width, 11. + off, 13. + off, 1.)

    tofDiff_e_mu = getTofDiff('e','mu',momentum)
    tofDiff_e_pi = getTofDiff('e','pi',momentum)
    tofDiff_e_p = getTofDiff('e','p',momentum)
    tofDiff_e_d = getTofDiff('e','d',momentum)

    print(f'ToF diffs for momentum {momentum}: mu-e: {tofDiff_e_mu:2.2f}, pi-e: {tofDiff_e_pi:2.2f}, p-e: {tofDiff_e_p:2.2f}, d-e: {tofDiff_e_d:2.2f}')

    if abs(momentum) < 330:
        print('Assuming low momemtum run, will try to fit e/mu/pi')
        fitmu = Fit(hTOFAll, '_mu', 11.7 + off + tofDiff_e_mu, width, 11. + off + tofDiff_e_mu, 13. + off + tofDiff_e_mu, 1.)
        fitpi = Fit(hTOFAll, '_pi', 11.7 + off + tofDiff_e_pi, width, 11. + off + tofDiff_e_pi, 13. + off + tofDiff_e_pi, 1.)
        stuff.append([fitmu, fitpi])

    elif abs(momentum) < 700:
        print('Assuming medium momemtum run, will try to fit e/mu+pi')
        fitmupi = Fit(hTOFAll, '_mupi', 11.7 + off + tofDiff_e_mu, width, 11. + off + tofDiff_e_mu, 13. + off + tofDiff_e_mu, 1.)
        stuff.append(fitmupi)
    else:
        print('Assuming low momemtum run, will try to fit e/p/d')
        fitp = Fit(hTOFAll, '_p',  16. + off, width, 14 + off,  17.5 + off, 1.)
        fitd = Fit(hTOFAll, '_d',  25. + off, 0.8, 24.3 + off, 25.7 + off, 0.8)
        stuff.append([fitp, fitd])
        
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


    
    #hTOFEl.SetLineColor(ROOT.kRed)
    #hTOFEl.SetLineWidth(2)
    #hTOFEl.Draw('hist same')

    ROOT.gPad.SetLogy(1)
    ROOT.gPad.Update()

    for can in cans:
        can.Print(pngdir + can.GetName() + '.png')
        can.Print(pdfdir + can.GetName() + '.pdf')
    
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

