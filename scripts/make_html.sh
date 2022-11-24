#!/bin/bash

mkdir TOFfit_Neg
mkdir TOFfit_Pos

cp output*n_*TOFfit*.png TOFfit_Neg/
cp output*p_*TOFfit*.png TOFfit_Pos/

for dd in TOFfit_Neg TOFfit_Pos ; do
cp scripts/CreateHTMLpage_noConv.sh ${dd}/
cd ${dd}
./CreateHTMLpage_noConv.sh
cd -
done



for i in `ls */*.html` ; do
  echo "google-chrome $i &"
done
