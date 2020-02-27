#!/bin/bash
#---------------------------------------------
# QCEIMS SLURM SCRIPT
#---------------------------------------------
#SBATCH --partition=production    # partition to submit to [production, intel]
###SBATCH --exclude=gaggle-[0,1]  # calling a specific gaggle node on [intel]
##SBATCH --job-name=work  # Job name
#SBATCH --nodes=1               # single node, anything more than 1 will not run
###SBATCH --exclusive             # mostly for use on intel partition
#SBATCH --ntasks=1         # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=2G      # memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-0:10:00     # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=plot%x_output.txt  # STDOUT
#SBATCH --error=plot%x_error.txt   # STDERR
#SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=sywang@ucdavis.edu
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
if [ ! -f /tmp/$user/qms$SLURM_JOB_NAME/qceims ]; then

	# create custom qceims for a specific user // folder name is hardcoded for each user
	mkdir -p /tmp/$user/qms$SLURM_JOB_NAME/

	# copy qceims work files to /qms
	cp -rT /share/fiehnlab/software/qceims4.0/ /tmp/$user/qms$SLURM_JOB_NAME/


	# change the executable mode
	chmod +x /tmp/$user/qms$SLURM_JOB_NAME/qceims

fi

cd /tmp/$user/qms$SLURM_JOB_NAME/
#set environment
source /share/fiehnlab/software/qceims4.0/.envrc

# make qceims, getres, plotms and mndo99 available
	export PATH=$PATH:/tmp/$user/qms$SLURM_JOB_NAME/


#pqceims results
DIR=TMPQCEIMS
cd $workdir$SLURM_JOB_NAME/
if [ ! -d "$DIR" ]; then
    echo "run qceims first!"
    exit
fi
cd $DIR

echo "start to read out/res file"
for vz in *
do
  cd $vz
  if [ -e  qceims.res ]
  then
  cat qceims.res >> ../../qceims.res
  cat qceims.out >> ../../qceims.out
  fi
  cd ..
done
cd $workdir$SLURM_JOB_NAME/
# get results - collect all
getres

# plotms - and print result.jdx with spectrum
plotms
sed -i "s/TITLE=/TITLE=${SLURM_JOB_NAME}"/ result.jdx
mv result.jdx GFn$SLURM_JOB_NAME.jdx
cp GFn$SLURM_JOB_NAME.jdx ../GFn$SLURM_JOB_NAME.jdx
mkdir arrayout
mv *-*.out arrayout
rm -rf /tmp/$user/qms$SLURM_JOB_NAME/
echo "Job ended at " `date`
