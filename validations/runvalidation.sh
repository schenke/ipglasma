#!/bin/bash
#SBATCH --job-name=validation
#SBATCH --account=<your_account>
#SBATCH --partition=small
#SBATCH --time=2:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000M

# Each event runs single-threaded; the workers run the events in parallel
export OMP_NUM_THREADS=1

# Place and bind threads to single cores
# Comment the following lines if binding is not desired
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
module add python-data gsl fftw

taskset -pc $$

# Run the program
python3 parallel_test_vector_meson_production.py --slurm \
    --max-workers "${SLURM_NTASKS:-1}" \
    --datadir /path/to/datadir \
    --maxevents 1500 \
    --subnucleondiffraction-path /path/to/subnucleondiffraction \
    --ipglasma-path /path/to/ipglasma

