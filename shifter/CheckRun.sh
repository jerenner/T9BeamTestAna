#!/bin/bash


if [ $# -lt 1 ] ; then
    echo "Usage $0 RUNNO"
    exit 1
fi



run=$1

echo "Thanks for checking end-of-run plots for run ${run} !"

infile=ntuple_files/ntuple_000${run}.root
if ! [ -f $infile ] ; then
    echo "ERROR: ntuple file $infile does not exist, maybe was not transferred from daq PC yet?"
    exit 1

fi


root -b -q -l "./macros/MakeDataPlots.C(\"${infile}\", 1000)"


hfile=histos/ntuple_000${run}_plots.root
if ! [ -f $hfile ] ; then
    echo "ERROR, it seems like the histogram file $hfile was not produce, did previous root macro succeed?"
    exit 1
fi

root -l "macros/quickPlots1d.C(\"${hfile}\")"
root -l "macros/morePlots.C(\"${hfile}\")"


