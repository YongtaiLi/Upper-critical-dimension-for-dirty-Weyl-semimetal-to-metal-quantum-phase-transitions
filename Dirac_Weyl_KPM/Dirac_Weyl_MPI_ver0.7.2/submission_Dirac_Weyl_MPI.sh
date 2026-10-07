#!/bin/bash

#SBATCH -A bir218_083127		# account (from which the SUs are granted)
#SBATCH -J High_Dim	 		# job name
#SBATCH -p hawkcpu			# queue (partition) that computes this job
#SBATCH -N 1				# total number of nodes, usually 1 (it is strongly discouraged to use > 2, even if your jobs are parallelized!!)
#SBATCH -n 1				# number of tasks that resources are allocated for
#SBATCH -c 1 				# number of processors (i.e., CPUs) per task
#SBATCH -t 48:00:00			# Max. Time consumption, capped at 72h. The shorter the time, the quicker a job gets scheduled to begin
#SBATCH -o datafiles/kwant-%j.out	# brief output and error file name ("%j" expands to jobID)
#SBATCH --mail-type=ALL 		# Send (ALL kinds of) update to the user's email
#SBATCH --mail-user=yol321@lehigh.edu	# the email of the user

###SBATCH -w hawk-a123 			# Optional. To specify the name of the node the job to submit to.

module load arch/cascade24v2		## module load arch/ice24v2, alternatively for some partitions with diff. architecture
module load miniconda3

conda activate ~/.conda/envs/kwant_environ  ## the location of conda environment that harbors python-kwant libraries. The name of the env is called "kwant_environ".

echo "Running on: $(hostname)"
echo "Python: $(which python)"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

## srun python Dirac_Weyl_main.py	## For MPI-parallelized job (Please use the serial-job command if SBATCH -n 1 )
python Dirac_Weyl_main.py		## For serial job

echo "    "
echo "Job done!"

echo "============== SLURM accounting info =============="

sacct -j "$SLURM_JOB_ID" -X \
        --format=JobID,JobName%16,Partition,Account%20,Elapsed,ElapsedRaw,AllocCPUS,AllocTRES%35

echo "    "
echo "---------------------------------------------------" ## Start calculating estimated SUs consumed.
echo "    "

sacct -j "$SLURM_JOB_ID" -X -n -P \
    --format=JobIDRaw,ElapsedRaw,AllocCPUS |
awk -F'|' -v job="$SLURM_JOB_ID" '
    $1 == job {
        printf "Approximate CPU service units: %.6f SU\n",
               $2 * $3 / 3600
    }
'

echo "==================================================="

exit

##########################################################
### To run, type:
### $ sbatch <name_of_this_submission_script>
##########################################################
