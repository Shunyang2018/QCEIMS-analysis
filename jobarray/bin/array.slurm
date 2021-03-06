#!/usr/bin/env bash
#---------------------------------------------
# QCEIMS Array SCRIPT
#---------------------------------------------
#SBATCH --partition=production     # partition to submit to [production, intel]
### SBATCH --nodelist=             # calling a specific gaggle node on [intel]
#SBATCH --job-name=om2-v12         # Job name
### SBATCH --nodes=1              # single node, anything more than 1 will not run
### SBATCH --exclusive             # mostly for use on intel partition
#SBATCH --ntasks=4                 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=4G                   # memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-00:20:00          # expected time of completion in hours, minutes, seconds, default 1-day
### SBATCH --output=Array.%A_%a.out  # STDOUT
### SBATCH --error=Array.%A_%a.error   # STDERR
### SBATCH --mail-user=you@gmail.com

workdir=$1
user=$2
# echo time and date
echo "Job started at " `date`
echo ''


# echo which node we are
echo  $HOSTNAME
echo ''

# echo the CPU
lscpu
echo ''

# echo the  job id
echo 'jobid: '${SLURM_JOB_ID}
echo ''


# check if qceims is already installed on the node, if not install it
if [ ! -f /tmp/$user/qms/qceims ]; then

	# create custom qceims for a specific user // folder name is hardcoded for each user
	mkdir -p /tmp/$user/qms/

	# copy qceims work files to /qms
	cp -rT /share/fiehnlab/software/qceims4.0/ /tmp/$user/qms


	# change the executable mode
	chmod +x /tmp/$user/qms/qceims

fi

cd /tmp/$user/qms
#set environment
ls -la
source /share/fiehnlab/software/qceims4.0/.envrc

# make qceims, getres, plotms and mndo99 available
	export PATH=$PATH:/tmp/$user/qms

# make new unique dir including username and  job-name
mkdir /tmp/$user/${SLURM_JOB_ID}/

# copy data from share to single node
cp -r $workdir$SLURM_JOB_NAME/TMPQCEIMS/TMP.${SLURM_ARRAY_TASK_ID} /tmp/$user/${SLURM_JOB_ID}/

# cd into same dir
cd /tmp/$user/${SLURM_JOB_ID}/TMP.${SLURM_ARRAY_TASK_ID}/

# run qceims production "> qceims.out 2>&1"
qceims -prod > qceims.out 2>&1

# wait
wait

# collect output back to share
cp -r /tmp/$user/${SLURM_JOB_ID}/TMP.${SLURM_ARRAY_TASK_ID}/* $workdir$SLURM_JOB_NAME/TMPQCEIMS/TMP.${SLURM_ARRAY_TASK_ID}

# remove the single directory
rm -rf /tmp/$user/${SLURM_JOB_ID}/TMP.${SLURM_ARRAY_TASK_ID}/
rm -rf /tmp/$user/qms
# tell its finished
# touch FINISHED

# echo that job ended
echo "Job ended at " `date`

