#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

#SBATCH --job-name=validation
#SBATCH --account=<your_account>
#SBATCH --partition=small
#SBATCH --time=2:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000M

# Each event runs single-threaded; the workers run the events in parallel.
# For n threads per event, set --cpus-per-task=n (and --ntasks to the number
# of events run at a time), OMP_NUM_THREADS=n, OMP_PROC_BIND=close instead of
# spread, and OMP_WAIT_POLICY=active: with --slurm, each event then runs as
# an srun step with n CPUs (see "Threads and MPI ranks" in the README).
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

