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
    }
    
    nChannels = 19 # 32
    Hs = []
    Txts = []
    ChNames =  {0: 'ACT-00', 1: 'ACT-01', 2: 'ACT-10', 3: 'ACT-11', 4: 'ACT-20', 5: 'ACT-21', 6: 'ACT-30', 7: 'ACT-31',
                8: 'TOF-00', 9: 'TOF-01', 10: 'TOF-10', 11: 'TOF-11', 12: 'TOF-20', 13: 'TOF-21', 14: 'TOF-30', 15: 'TOF-31',
                16: 'HC-00', 17: 'HC-10', 18: 'LG-00',
                19: 'X', 20: 'X', 21: 'X', 22: 'X', 23: 'X',
                24: 'X', 25: 'X', 26: 'X', 27: 'X', 28: 'X', 29: 'X', 30: 'X', 31: 'X' }


    ftag = filename.split('/')[-1].replace('output_','').replace('_plots.root','')

    os.system('mkdir -p pdf png')
    
    for hbasename in hbasenames:

        hs = []
        txts = []


        canname = 'WCTEJuly2023_Quick1D_{}_{}'.format(ftag, hbasename)
        canname = canname.replace('_list_root','')

        can = ROOT.TCanvas(canname, canname, 0, 0, 1600, 800)
        cans.append(can)
        can.Divide(8,4)
        #can.Divide(4,2)
        
        for ich in range(0, nChannels):
            # hack just for old tofs
            #if not ( ich >= 8 and ich <= 15):
            #    continue
            hname = hbasename + str(ich)
            h = rfile.Get(hname)
            print('Pushing ', ich, hname)
            hs.append(h)
        Hs.append(hs)
        
        for h in hs:
            ih = hs.index(h)
            can.cd(ih+1)
            h.SetStats(0)
            #if not 'Time' in h.GetName():
            #    ROOT.gPad.SetLogy(1)
            #h.GetYaxis().SetRangeUser(1.e-4, h.GetYaxis().GetXmax())
            h.SetFillColor(hbasenames[hbasename])
            h.SetFillStyle(1111)
            h.SetTitle(ChNames[hs.index(h)])

            gmeans = []
            gsigmas = []
            xmaxes = []
            
            if 'Time' in h.GetName() and ih >= 8 and ih <= 15:
                #print('*** fitting {}'.format(h.GetName()))
                hname = h.GetName()
                fname = 'fit{}'.format(ih)
                fit = ROOT.TF1(fname, '[0]*exp(-(x-[1])^2/(2*[2]^2))', 0., 520.)
                fit.SetParameters(h.GetMaximum()/6., 80., 5.)
                h.Fit(fname, 'q', '')
                mean = fit.GetParameter(1)
                sigma = fit.GetParameter(2)
                #print(f'1) {hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
                sf = 2.
                h.Fit(fname, 'q', '', mean - sf*sigma, mean + sf*sigma)
                mean = fit.GetParameter(1)
                sigma = fit.GetParameter(2)
                print(f'{hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
                h.SetMaximum(1.2*h.GetMaximum())
                h.Draw('hist')
                h.GetXaxis().SetRangeUser(50, 120)
                fit.Draw('same')
                chi2 = fit.GetChisquare()
                ndf = fit.GetNDF()
                maxb = h.GetMaximumBin()
                bw = h.GetBinWidth(maxb)
                ymax = h.GetBinContent(maxb)
                xmax = h.GetBinCenter(maxb)
                if ndf > 0:
                    print(f'  chi2/ndf={chi2/ndf:1.3f}, bw={bw} mean-ymax={mean-xmax:1.3f}')
                stuff.append(fit)
                gmeans.append(mean)
                gsigmas.append(sigma)
                xmaxes.append(xmax)
            else:
                h.Draw('hist')

            print(gmeans)
            print(gsigmas)
            print(xmaxes)

            #mean
            
            # todo: times subtractions 
            
            #txt= 'Ch {} {} p={} MeV/c'.format(hs.index(h)+1, basetag, ftag.replace('n','').replace('p',''))
            #if 'p' in ftag:
            #    txt = txt + ' (Pos)'
            #elif 'n' in ftag:
            #    txt = txt + ' (Neg)'
            #else:
            #    txt = txt + ' (undef)'
            #text = ROOT.TLatex(0.120, 0.93, txt)
            #text.SetNDC()
            #text.Draw()
            #txts.append(text)


    # 2023 toch channels time diff analysis / calibration
    htdiffnames = []
    base = 'hTimeDiffTOF'
    htofdiffs = {}
    tofmeans = {}
    for itof in range(0,2):
        htofdiffs[itof] = []
        tofmeans[itof] = []
        for itch in range(1,4):
            hname = base + str(itof) + str(itch)
            htdiffnames.append(hname)
            h = rfile.Get(hname)
            htofdiffs[itof].append(h)

            
    canname = 'TofDiffsWCTEJuly2023_Quick1D_{}_{}'.format(ftag, hbasename)
    canname = canname.replace('_list_root','')
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    can.Divide(3,2)
    cans.append(can)
    ic = 1

    for itof,hs in htofdiffs.items():
        ih = -1
        for h in hs:
            ih = ih + 1
            can.cd(ic)
            h.Draw('hist')
            #h.GetXaxis().SetRangeUser(-4,8.)
            fname = 'fit_tofs_{}'.format(ic)
            fit = ROOT.TF1(fname, '[0]*exp(-(x-[1])^2/(2*[2]^2))', h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
            fit.SetParameters(h.GetMaximum()/6., h.GetMean(), h.GetStdDev())
            h.Fit(fname, 'q', '', )
            fit.Draw('same')
            mean = fit.GetParameter(1)
            sigma = fit.GetParameter(2)
            print(f'tof fit {itof} {ih}: {hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
            stuff.append(fit)
            ic = ic+1
            tofmeans[itof].append(mean)
    print(tofmeans)
    for can in cans:
        can.Update()
        can.Print('png/' + can.GetName() + '.png')
        can.Print('pdf/' + can.GetName() + '.pdf')
     
   # if not gBatch:
   #     ROOT.gApplication.Run()
   # return


    ######## absolute TOF measurement ##########
    htTOFnames = []
    base = 'hTimeTOF'
    htofTOFs = {}
    htofTOFs[0] = []
    for itch in range(0,4):
        hname = base + str(itch)
        htTOFnames.append(hname)
        h = rfile.Get(hname)
        htofTOFs[0].append(h)
    
    canname = 'AbsoluteTof_{}'.format(ftag[10:])
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    can.Divide(2,2)
    cans.append(can)
    ic = 1

    for itof, hs in htofTOFs.items():
        ih = -1
        for h in hs:
            ih = ih + 1
            can.cd(ic)
            ROOT.gPad.SetLogy(1);
            print("ok")
            h.Draw('hist')
            
            #h.GetXaxis().SetRangeUser(-4,8.)
            fname = 'fit_tofs_{}'.format(ic)
            fit = ROOT.TF1(fname, '[0]*exp(-(x-[1])^2/(2*[2]^2)) + [3]*exp(-(x-[4])^2/(2*[5]^2))', h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
            fit.SetParameters(h.GetMaximum()/6., h.GetMean()-1, h.GetStdDev(), h.GetMaximum()/6., h.GetMean()+1, h.GetStdDev())
            #try to fiut the peaks one by one
            #g1 = ROOT.TF1("g1", "gaus", 0, 15);
            #g2 = ROOT.TF1("g2", "gaus", 12, 20);
            #g3 = ROOT.TF1("g3", "gaus", 25, 30);

            #h.Fit(g1, "R");
            #h.Fit(g2, "R+");
            #h.Fit(g3, "+", "", 0, 50);


            #fit = g1 + g2 + g3 
            h.Fit(fname, 'q', '', )
            fit.Draw('same')
            #g1.Draw('same')
            #g2.Draw('same')
            #g3.Draw('same')
            mean = fit.GetParameter(1)
            sigma = fit.GetParameter(2)
            mean2 = fit.GetParameter(4)
            sigma2 = fit.GetParameter(5)
            mean3 = fit.GetParameter(6)
            sigma3 = fit.GetParameter(7)
            print("ok")
            print(f'tof fit {itof} {ih}: {hname } mean={mean:1.3f} ns; sigma={sigma:1.3f} ns')
            stuff.append(fit)
            ic = ic+1
            tofmeans[itof].append(mean)

    
##################################s
#   plot 2d hist act vs lead     #
##################################

    canname = 'ACT2+3vsLeadGlass_{}'.format(ftag[10:])
    hname = 'hRef_pbA_act23A'
    h = rfile.Get(hname)
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    cans.append(can)
    can.cd()
    integral_full = h.Integral()
    h.GetXaxis().SetRangeUser(0, 8)
    h.GetYaxis().SetRangeUser(0, 30)
    x1, x2, y1, y2 = 0.2, 1.3, 1.5, 4.5
    bx1 = h.GetXaxis().FindBin(x1)
    bx2 = h.GetXaxis().FindBin(x2)
    by1 = h.GetYaxis().FindBin(y1)
    by2 = h.GetYaxis().FindBin(y2)
    integral_zoom = h.Integral(bx1, bx2, by1, by2)
    print(f'integral: ful={integral_full}, zoom={integral_zoom}')
    h.SetTitle('ACT2+3 vs Lead Glass {}'.format(ftag[10:]))
    h.Draw("colz")
    adjustStats(h)
    
    if False:
        line = makeLine(x1, x2, y1, y1)
        lines.append(line)
        line = makeLine(x1, x2, y2, y2)
        lines.append(line)
        line = makeLine(x1, x1, y1, y2)
        lines.append(line)
        line = makeLine(x2, x2, y1, y2)
        lines.append(line)
        txt = ROOT.TLatex(0.14, 0.8, f'box: {integral_zoom/1000:1.2f}k evts, full: {integral_full/1000:1.2f}k evts')
        txt.SetTextSize(0.04)
        txt.SetNDC()
        txt.Draw()
        stuff.append(txt)


    canname = 'TOFvsACT3_{}'.format(ftag[10:])
    hname = 'hRef_TOFACT3A'
    h = rfile.Get(hname)
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    cans.append(can)
    can.cd()
    h.SetTitle('TOF vs ACT3 {}'.format(ftag[10:]))
    h.GetXaxis().SetRangeUser(0, 4)
    h.GetYaxis().SetRangeUser(0, 30)
    h.Draw("colz")
    adjustStats(h)

    #h.colz()

##################################
#       plots all the canvas     #
##################################

    print(tofmeans)

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

