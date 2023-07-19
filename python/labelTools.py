#!/usr/bin/python

import ROOT

from data_runs import *

ChNames =  {0: 'ACT-00', 1: 'ACT-01', 2: 'ACT-10', 3: 'ACT-11', 4: 'ACT-20', 5: 'ACT-21', 6: 'ACT-30', 7: 'ACT-31',
            8: 'TOF-00', 9: 'TOF-01', 10: 'TOF-10', 11: 'TOF-11', 12: 'TOF-20', 13: 'TOF-21', 14: 'TOF-30', 15: 'TOF-31',
            16: 'HC-00', 17: 'HC-10', 18: 'LG-00',
            19: 'X', 20: 'X', 21: 'X', 22: 'X', 23: 'X',
            24: 'X', 25: 'X', 26: 'X', 27: 'X', 28: 'X', 29: 'X', 30: 'X', 31: 'X' }


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
    
