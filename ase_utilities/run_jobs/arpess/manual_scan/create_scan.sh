#!/bin/bash

#-0.331 0.019 0.05


for i in $(seq $1 $3 $2)
  do echo "creating scan_$i"
  mkdir scan_$i
  cd scan_$i
  cp ../{in.traj,run.py,vasp.sh} .
  sed -i "s/initial_value=0.0/initial_value=$i/g" run.py
  new_submission vasp.sh
  cd ..
  done
