#!/usr/bin/env bash
#---------------------------------------------
# QCEIMS Array SCRIPT
#---------------------------------------------
#SBATCH --partition=intel   # partition to submit to [production, intel]
### SBATCH --nodelist=             # calling a specific gaggle node on [intel]
###SBATCH --exclude=gaggle-[0,1]
#SBATCH --job-name=om2-v12         # Job name
#SBATCH --nodes=1              # single node, anything more than 1 will not run
#SBATCH --ntasks=1                 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=8G                   # memory pool all cores, default is 2GB per cpu
#SBATCH --time=7-960:00:00          # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=Array.%A_%a.out  # STDOUT
#SBATCH --error=Array.%A_%a.error   # STDERR
##SBATCH --mail-user=sywang@ucdavis.edu


unlog
aklog
module load orca/4.0.0.2
rm *.xyz
rm qceims.out
rm qceims.res
rm ready
rm trj.*
echo 'orca' > qceims.in
echo 'ip-orca' >> qceims.in
echo 'pbe0' >> qceims.in
echo 'SV(P)' >> qceims.in
/share/fiehnlab/software/qceims4.0/qceims -p > qceims.out

