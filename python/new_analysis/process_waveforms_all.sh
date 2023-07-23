#!/bin/bash
shopt -s extglob
datadir=/media/wcte/T7/WCTE/data/2023/
dir=root_files/
for i in `cd $datadir/$dir ; ls root_run_*.root` ; do
  outfile="${datadir}/ntuple_files/${i/root_run/ntuple}"
#  if ! [ -f $outfile ] ; then
    script="python3 /home/wcte/np/T9BeamTestAna/python/new_analysis/process_waveform_analysis.py"
    config="/home/wcte/np/T9BeamTestAna/config/config.json"
    echo "${script} ${datadir}/${dir}/${i} ${config} ${outfile}"
    $script "${datadir}/${dir}/${i}" "${config}" "${outfile}"
#  else
#    echo "NOT processing as output ${outfile} exist!"
#  fi
done
