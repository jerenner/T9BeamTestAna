#!/snap/bin/pyroot
# was: #!/usr/bin/python3
# Pá 14. července 2023, 18:39:38 CEST

#from __future__ import print_function

from data_runs import *
from tofUtil import *
from labelTools import *

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
def makeLines(h, eoff, parts, momentum):
    lines = []
    #y1 = h.GetYaxis().GetXmin()
    #y2 = h.GetYaxis().GetXmax()
    y1 = h.GetMaximum()
    y2 = h.GetMinimum()
    te = getTof(ms['e'], momentum) + eoff
    for part in parts:
        dt = getTofDiff('e', part, momentum)
        print(f'makeLines {part}: dt={dt} ns')
        print('line coors: ', te + dt, y1, te + dt, y2)
        line = ROOT.TLine(te + dt, y1, te + dt, y2)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.SetLineColor(pcols[part])
        line.Draw()
        lines.append(line)
    return lines


##########################################
# cts     ...  central peak time of assumed gauss
# w       ... width
# t1, t2: ... fit window

def Fit(h, tag, momentum, ct, w, t1, t2, peaksf = 1.):
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
    h.Fit(fname, '', '', t1, t2)
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

    te = getTof(ms['e'], momentum)
    eoff = fit.GetParameter(1) - te
    
    parts = ['p', 'd']
    lines = makeLines(h, eoff, parts, momentum)
    stuff.append(lines)    
    
    return fit


##########################################
# cts    ... central peak time of assuemd gauss
# w      ... widths
# t1, t2 ... full fit window

def MultiFit(h, tag, momentum, cts, ws, t1, t2, peaksfs = [1., 1., 1.]):
    fname = 'fit{}'.format(tag)
    hname = h.GetName()
    fitform = ''
    ngpars = 3
    # possible future issue: same sigma of pi and mu?
    # fragile change to make...
    
    for ifit in range(0,len(cts)):
        fitform = fitform + '[{}]*exp(-(x-[{}])^2/(2*[{}]^2))'.format(ifit*ngpars, ifit*ngpars + 1, ifit*ngpars + 2 )
        if ifit < len(cts)-1:
            fitform = fitform + ' + '
    print('Fit formula: ', fitform)
    fit = ROOT.TF1(fname, fitform, t1, t2)
    ampl = h.GetBinContent(h.FindBin(cts[0]))
    for ifit in range(0,len(cts)):
        fit.SetParameter(ifit*ngpars+1, ampl / peaksfs[ifit])
        fit.SetParameter(ifit*ngpars+1, cts[ifit])
        fit.SetParameter(ifit*ngpars+2, ws[ifit])
    for ipar in range(0,fit.GetNpar()):
        print('par{} initially {:1.2f}'.format(ipar, fit.GetParameter(ipar)))
    print('*** Fitting!')
    h.Fit(fname, '', '' , t1, t2)
    fit.SetLineColor(ROOT.kBlack)
    fit.SetLineStyle(2)
    #fit.SetLineColor(h.GetLineColor())
    fit.SetLineStyle(2)
    fit.Draw('same')

    eoff = fit.GetParameter(1)

    fits = [fit]
    for ifit in range(0,len(cts)):
        sfname = fname + str(ifit)
        fitform = '[{}]*exp(-(x-[{}])^2/(2*[{}]^2))'.format(0, 1, 2)
        single_fit = ROOT.TF1(sfname, fitform, t1, t2)
        for ip in range(0, ngpars):
            single_fit.SetParameter(ip, fit.GetParameter(ifit*ngpars + ip) )
        fits.append(single_fit)
    for sfit in fits[1:]:
        sfit.SetLineStyle(2)
        sfit.SetLineColor(ROOT.kBlue)
        sfit.Draw('same')


    te = getTof(ms['e'], momentum)
    eoff = fit.GetParameter(1) - te

    parts = ['e', 'mu', 'pi']
    lines = makeLines(h, eoff, parts, momentum)
    stuff.append(lines)
    
    return fits




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
    srun = ''
    tokens = fileName.split('/')[-1].split('_')
    for token in tokens:
        if '00' in token:
            srun = token.replace('000','')
        
    momentum = getMomentum(srun)
    print(f'Assuming run {srun} and momentum {momentum}')

    suff = ''
    if abs(momentum) < 500:
        suff = 'Low'
        print('Low momentum run, looking at zoomed version of tof histos!')
    hTOFAll = inFile.Get("hTOFOther" + suff)
    hTOFEl = inFile.Get("hTOFEl" + suff)


    signedmomentum = str(abs(momentum))
    if momentum > 0:
        signedmomentum = signedmomentum + 'Pos'
    else:
        signedmomentum = signedmomentum + 'Neg'
    canname = 'FitToF_run{}_{}'.format(srun, signedmomentum)
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

    tofDiff_e_mu = getTofDiff('e','mu',momentum)
    tofDiff_e_pi = getTofDiff('e','pi',momentum)
    tofDiff_e_p = getTofDiff('e','p',momentum)
    tofDiff_e_d = getTofDiff('e','d',momentum)

    print(f'ToF diffs for momentum {momentum}: mu-e: {tofDiff_e_mu:2.2f}, pi-e: {tofDiff_e_pi:2.2f}, p-e: {tofDiff_e_p:2.2f}, d-e: {tofDiff_e_d:2.2f}')

    if abs(momentum) < 300:
        print('Assuming low momentum run, will try to fit e/mu/pi')
        tcs = [11.7, 11.7 + off + tofDiff_e_mu, 11.7 + off + tofDiff_e_pi]
        ws = [0.25, 0.3, 0.35]
        sfs = [1., 10., 25.]
        fits = MultiFit(hTOFAll, '_mupi', momentum, tcs, ws, 11., 14., sfs)
        stuff.append(fits)
    elif abs(momentum) < 700:
        print('Assuming medium momentum run, will try to fit e/mu+pi')
        tcs = [11.7, 11.7 + off + tofDiff_e_mu]
        ws = [0.25, 0.3]
        sfs = [1., 5.]
        fits = MultiFit(hTOFAll, '_mupi', momentum, tcs, ws, 11., 14., sfs)
        stuff.append(fits)
    else:
        print('Assuming high momentum run, will try to fit e/p/d')
        fite = Fit(hTOFAll, '_el', momentum, 11.7 + off, width, 11. + off, 13. + off, 1.)
        
        fitp = Fit(hTOFAll, '_p', momentum,  16. + off, width, 14 + off,  17.5 + off, 1.)
        fitd = Fit(hTOFAll, '_d', momentum,  25. + off, 0.8, 24.3 + off, 25.7 + off, 0.8)
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

    ROOT.gPad.Update()

    pnote = makeMomentumLabel(srun, 0.14, 0.92)
    stuff.append(pnote)
    for can in cans:
        can.cd()
        pnote.Draw()
        can.Print(pngdir + can.GetName() + '_liny.png')
        can.Print(pdfdir + can.GetName() + '_liny.pdf')
        ROOT.gPad.SetLogy(1)
        can.Print(pngdir + can.GetName() + '_logy.png')
        can.Print(pdfdir + can.GetName() + '_logy.pdf')

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

