#!/usr/bin/python

import ROOT

from data_runs import *

def makeMomentumLabel(srun, x = 0.14, y = 0.84, size = 0.04):
    momentum = getMomentum(srun)
    tag = 'Pos'
    if momentum < 0:
        tag = 'Neg'
    pnote = ROOT.TLatex(x, y, f'WCTE TB2023 run {srun}, p={abs(momentum)} MeV/c {tag}')
    pnote.SetTextSize(size)
    pnote.SetNDC()
    return pnote


def adjustStats(h):
    ROOT.gPad.Update()
    st = h.GetListOfFunctions().FindObject("stats")
    st = ROOT.gPad.GetPrimitive("stats")
    st.SetX1NDC(0.7)
    st.SetX2NDC(0.9)
    st.SetY1NDC(0.65)
    st.SetY2NDC(0.9)
    
