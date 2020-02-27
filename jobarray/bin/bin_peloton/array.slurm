#!/usr/bin/env bash
#---------------------------------------------
# QCEIMS Array SCRIPT
#---------------------------------------------
#SBATCH --partition=med     # partition to submit to [production, intel]
#SBATCH --ntasks=8                 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=20G                   # memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-4:20:00          # expected time of completion in hours, minutes, seconds, default 1-day
# SBATCH --output=Array.%A_%a.out  # STDOUT
# SBATCH --error=Array.%A_%a.error   # STDERR

workdir=$1
user=$2
in=$3
# echo time and date
echo "Job started at " `date`
echo ''
export OMP_NUM_THREADS=8,1

# echo which node we are
echo  $HOSTNAME
echo ''
module load orca/3.0.3

# echo the  job id
echo 'jobid: '${SLURM_JOB_ID}
echo ''


# check if qceims is already installed on the node, if not install it
if [ ! -f /tmp/$user/qms$SLURM_JOB_NAME${SLURM_ARRAY_TASK_ID}/qceims ]; then

	# create custom qceims for a specific user // folder name is hardcoded for each user
	mkdir -p /tmp/$user/qms$SLURM_JOB_NAME${SLURM_ARRAY_TASK_ID}/

	# copy qceims work files to /qms
	cp -rT /home/shunyang/qceims4.0/ /tmp/$user/qms$SLURM_JOB_NAME${SLURM_ARRAY_TASK_ID}/


	# change the executable mode
	chmod +x /tmp/$user/qms$SLURM_JOB_NAME${SLURM_ARRAY_TASK_ID}/qceims

fi

cd /tmp/$user/qms$SLURM_JOB_NAME${SLURM_ARRAY_TASK_ID}/
#set environment
ls -la
source /home/shunyang/qceims4.0/.envrc

# make qceims, getres, plotms and mndo99 available
export PATH=$PATH:/tmp/$user/qms$SLURM_JOB_NAME${SLURM_ARRAY_TASK_ID}

# make new unique dir including username and  job-name
mkdir /tmp/$user/${SLURM_JOB_ID}/

# copy data from share to single node
cp -r $workdir$SLURM_JOB_NAME/TMPQCEIMS/TMP.${SLURM_ARRAY_TASK_ID} /tmp/$user/${SLURM_JOB_ID}/

# cd into same dir
cd /tmp/$user/${SLURM_JOB_ID}/TMP.${SLURM_ARRAY_TASK_ID}/
echo -e $in > qceims.in 
# run qceims production "> qceims.out 2>&1"
qceims -prod > qceims.out 2>&1

# wait
wait

# collect output back to share
cp -r /tmp/$user/${SLURM_JOB_ID}/TMP.${SLURM_ARRAY_TASK_ID}/* $workdir$SLURM_JOB_NAME/TMPQCEIMS/TMP.${SLURM_ARRAY_TASK_ID}

# remove the single directory
rm -rf /tmp/$user/qms$SLURM_JOB_NAME${SLURM_ARRAY_TASK_ID}/
rm -rf /tmp/$user/${SLURM_JOB_ID}/TMP.${SLURM_ARRAY_TASK_ID}/
find /tmp/$user/ -type d -empty -delete
# tell its finished
# touch FINISHED

# echo that job ended
echo "Job ended at " `date`

