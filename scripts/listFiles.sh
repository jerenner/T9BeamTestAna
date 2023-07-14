#!/bin/bash

mkdir -p lists
path=data/

for i in `cd $path ; ls root_run_*.root` ; do
    j=`basename $i .root`
    echo $j
    ls ${path}/${i} > lists/list_${j}.txt
done


