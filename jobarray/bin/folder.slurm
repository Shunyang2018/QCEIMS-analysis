#!/bin/bash
#---------------------------------------------
# QCEIMS SLURM SCRIPT
#---------------------------------------------
#SBATCH --partition=production    # partition to submit to [production, intel]
###SBATCH --exclude=gaggle-[0,1]  # calling a specific gaggle node on [intel]
##SBATCH --job-name=work  # Job name
#SBATCH --nodes=1               # single node, anything more than 1 will not run
###SBATCH --exclusive             # mostly for use on intel partition
#SBATCH --ntasks=16         # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=120G      # memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-2:30:00     # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=%x_output.txt  # STDOUT
#SBATCH --error=%x_error.txt   # STDERR
#SBATCH --mail-user=sywang@ucdavis.edu

#---------------------------------------------
# Example script for running QCEIMS on cluster
#---------------------------------------------
workdir=$1
user=$2
# echo time and date
echo "Job started at " `date`

# show hostname
echo $HOSTNAME

# show CPU info
lscpu

# create custom qceims for a specific user // folder name is hardcoded for each user
mkdir -p /tmp/$user/qms/

# copy qceims work files to /qms
# cp: omitting directory /share/fiehnlab/software/qceims/cluster/original-input
cp -rT /share/fiehnlab/software/qceims4.0/ /tmp/$user/qms

# change the executable mode
# ln -s /tmp/$user/qms/mndo99_20131011_intel64_composerxe-2013.4.183 /tmp/$user/qms/mndo99
# chmod +x /tmp/$user/qms/mndo99
chmod +x /tmp/$user/qms/qceims

# make qceims, getres, plotms and mndo99 available
export PATH=$PATH:/tmp/$user/qms

cd /tmp/$user/qms
#set environment
source /share/fiehnlab/software/qceims4.0/.envrc
# make directory for workfiles
mkdir -p /tmp/$user/$SLURM_JOB_NAME

# copy inchikey
# cp /share/fiehnlab/users/shunyang/qceims/project/structure/work-inchikey.txt /tmp/$user/work/

# copy mol
cp $workdir$SLURM_JOB_NAME/structure.tmol /tmp/$user/$SLURM_JOB_NAME

# change to workdir on tmp
cd /tmp/$user/$SLURM_JOB_NAME

# copy tmol to file "coord"
cp *.tmol coord

# start qceims MD (single core)
qceims

# run DFT in ORCA (single core)
qceims

# pack the TMPQCEIMS folder and remove it
cd /tmp/$user/
tar -zcvf $SLURM_JOB_NAME.tar.gz $SLURM_JOB_NAME

# move result file back to share
mv $SLURM_JOB_NAME.tar.gz $workdir

# clean up old qms dir and work folder
rm -rf /tmp/$user/qms
rm -rf /tmp/$user/$SLURM_JOB_NAME

#back to share
cd $workdir
tar xvzf $SLURM_JOB_NAME.tar.gz
rm $SLURM_JOB_NAME.tar.gz


# echo that job ended
echo "Job ended at " `date`


