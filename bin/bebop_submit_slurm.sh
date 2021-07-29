#!/bin/bash
#SBATCH --job-name=quopt7
#SBATCH --account=qc_err_sim
#SBATCH --partition=bdwall
##SBATCH --partition=bdws
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --output=quopt.out
#SBATCH --error=quopt.error
#SBATCH --mail-user=mbowman@anl.gov # Optional if you require email
#SBATCH --mail-type=ALL # Optional if you require email
#SBATCH --time=24:00:00
# Set up Environment
export PATH=/soft/ceres-solver/2.0.0/conda/bin:$PATH
# Test
echo -e "Slurm job ID: $SLURM_JOBID"
# cd $PBS_O_WORKDIR
cd $SLURM_SUBMIT_DIR
# A little useful information for the log file...
echo -e "Master process running on: $HOSTNAME"
echo -e "Directory is:  $PWD"
# Put in a timestamp
echo Starting executation at: `date`
cmd="srun mpirun -np 36 ./main"
echo The command is: $cmd
echo End PBS script information.
echo -e "All further output is from the process being run and not the pbs script.\n$cmd\n\n"
export OMP_NUM_THREADS=10
ls
pwd
$cmd
# Print the date again -- when finished
echo Finished at: `date`
