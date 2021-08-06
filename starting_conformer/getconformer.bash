#!/usr/bin/env bash
#author SYWANG
mkdir startconformer
DIR=/Users/shunyang/project/TMS/TMS/TMS-traj/4-chlorobutan-1-ol/TMPQCEIMS
if [ ! -d "$DIR" ]; then
    echo "run qceims first!"
    exit
fi

cd $DIR

for vz in TMP.*
do

  cd $vz
  pwd
  vz=${vz#*.}
  n=$(head -1 trj.$vz.1)

  echo $(($n+2))
  head -n $(($n+2)) trj.$vz.1 > ../../startconformer/trj.$vz.xyz
  cd ..
done

echo 'normal termination'
