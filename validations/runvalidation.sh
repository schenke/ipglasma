#!/bin/bash
#SBATCH --job-name=validation
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=2:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000M


# Set the number of threads based on cpus-per-task
export OMP_NUM_THREADS=1
#${SLURM_CPUS_PER_TASK:-1}


# Place and bind threads to single cores
# Comment the following lines if binding is not desired
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
module add python-data gsl fftw

taskset -pc $$

# Run the program
python3 parallel_test_vector_meson_production.py -max_workers "${SLURM_NTASKS:-1}" -datadir /path/todatadir -maxevents 1500 -subnucleondiffraction_path /projappl/lappi/heikki/subnucleondiffraction -ipglasma_path /scratch/lappi/heikki/ipglasma_upstream

