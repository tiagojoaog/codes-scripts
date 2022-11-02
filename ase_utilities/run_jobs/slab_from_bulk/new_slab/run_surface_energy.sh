#!/bin/bash

trajectory=$(ls slab_*traj)

for i in $trajectory;
do folder_name=$(echo $i | sed "s/.traj//g")
mkdir $folder_name
cd $folder_name
cp ../{run.py,vasp.sh} .
mv ../$i in.traj
sbatch vasp.sh
cd ..
done
