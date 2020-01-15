#!/usr/bin/env bash
#author SYWANG
###reference: https://hpc.nih.gov/docs/job_dependencies.html
###reference: https://stackoverflow.com/questions/52248393/having-a-job-depend-on-an-array-job-in-slurm

#use this code by 'bash array.bash STUCTURE.tmol'
file=$1
if [ ! $file ]; then
        echo "Usage: bash array.bash STUCTURE.tmol"
        exit
fi

#########config part
#put this file and bin folder and your structure *.tmol under workdir
workdir=/share/fiehnlab/users/shunyang/qceims/array/
#the name used in /tmp/user/qme
user=swang
#i didn't have a config about these: you can change it in the slurm file
#task       ncpu    mem
#folder     16      120G
#each array 4       32G
#plot        4       8G
#########

echo $file
echo workdir: $workdir
echo user:$user
path=`echo $file|sed 's/.tmol//g'`
mkdir $path
anumber=`grep -v '^\s*$' $file | wc -l`
#get number of trajectories by default setting, anumber*25
trajnumber=$(($(($anumber-2))*25))

echo 'array_number'
echo $trajnumber
#collect slurm files and structures
cp ${workdir}bin/folder.slurm $workdir$path/folder.slurm
cp ${workdir}bin/array.slurm $workdir$path/array.slurm
cp ${workdir}bin/plot.slurm $workdir$path/plot.slurm
cp $workdir$1 $workdir$path/structure.tmol
#go to the subfolder
cd $path
echo 'workpath'
pwd


echo 'sub job1'
jid1=$(sbatch --job-name=$path --parsable folder.slurm $workdir $user)
echo $jid1
echo 'sub job2'
jid2=$(sbatch  --dependency=afterany:$jid1 --job-name=$path --array=1-$trajnumber --parsable array.slurm $workdir $user)
echo $jid2
echo 'sub job3'
sbatch --dependency=singleton --job-name=$path plot.slurm $workdir $user

