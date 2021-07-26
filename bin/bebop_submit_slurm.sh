#!/bin/bash

#SBATCH --job-name=quopt
#SBATCH --account=qc_err_sim
##SBATCH --partition=bdwall
#SBATCH --partition=bdws
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=quopt.out
#SBATCH --error=quopt.error
#SBATCH --mail-user=mbowman@anl.gov # Optional if you require email
#SBATCH --mail-type=ALL # Optional if you require email
#SBATCH --time=1:00:00

# Set up Environment
export PATH=/soft/ceres-solver/2.0.0/conda/bin:$PATH

#Test
echo -e "Slurm job ID: $SLURM_JOBID"

#cd $PBS_O_WORKDIR
cd $SLURM_SUBMIT_DIR

# A little useful information for the log file...
echo -e "Master process running on: $HOSTNAME"
echo -e "Directory is:  $PWD"

# Put in a timestamp
echo Starting executation at: `date`

srun hostname | sort -u > node_list

#Add manager node to machinefile
head -n 1 node_list > machinefile.$SLURM_JOBID
#Add worker node to machinefile
cat node_list >> machinefile.$SLURM_JOBID

cmd="srun mpirun -np 1 ./main"

echo The command is: $cmd
echo End PBS script information.
echo -e "All further output is from the process being run and not the pbs script.\n$cmd\n\n"

export OMP_NUM_THREADS=10
export SLURM_HOSTFILE=machinefile.$SLURM_JOBID
ls
pwd
$cmd

# Print the date again -- when finished
echo Finished at: `date`
