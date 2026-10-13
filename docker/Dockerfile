# Image with IP-Glasma built from a branch of the GitHub repository, with
# the files to run it with Singularity/Apptainer.
#
# Build (from the repository root):
#   docker build -t ipglasma docker/
#   docker build -t ipglasma --build-arg IPGLASMA_BRANCH=devel docker/
# Run (the input's data files, e.g. nucleusQsTableFileName and
# nuclearConfigurationsPath, are relative to the working directory; the
# repository with the nucleus configurations is in /home/ipglasma):
#   docker run --rm -it -v "$PWD:/work" -w /work ipglasma ipglasma input
FROM ubuntu:24.04

LABEL maintainer="Bjoern Schenke <bschenke@bnl.gov>"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        git \
        libfftw3-dev \
        libgsl-dev \
        libopenmpi-dev \
        openmpi-bin \
        python3 \
        python3-h5py \
        python3-matplotlib \
        python3-numpy \
        python3-scipy \
        time \
        vim \
        wget \
        && \
    rm -rf /var/lib/apt/lists/*

ENV HOME=/home
ENV PATH="/home/ipglasma:${PATH}"
# Open MPI refuses to run as root, the default user in Docker
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# download and compile IP-Glasma and download the nucleus configurations;
# the executable ipglasma is installed in /home/ipglasma, and the source
# stays there since ipglasma calls utilities/combine_events_into_hdf5.py
# for writeOutputsToHDF5
ARG IPGLASMA_BRANCH=main
# number of parallel compile jobs; each one can take about 1 GB of memory
ARG BUILD_JOBS=4
RUN cd ${HOME} && \
    git clone --depth=1 https://github.com/schenke/ipglasma -b ${IPGLASMA_BRANCH} && \
    cd ipglasma && \
    cmake -B build && \
    cmake --build build -j "${BUILD_JOBS}" && \
    cmake --install build && \
    rm -rf build && \
    cd nucleusConfigurations && \
    bash download_nucleusTables.sh

# mount points of the Open Science Grid
RUN for MNTPOINT in \
        /cvmfs \
        /hadoop \
        /hdfs \
        /lizard \
        /mnt/hadoop \
        /mnt/hdfs \
        /xenon \
        /scratch \
        /spt \
        /stash2 \
    ; do \
        mkdir -p $MNTPOINT ; \
    done

# Singularity/Apptainer files
COPY .singularity.d /.singularity.d
RUN cd / && \
    ln -s .singularity.d/actions/exec .exec && \
    ln -s .singularity.d/actions/run .run && \
    ln -s .singularity.d/actions/shell .shell && \
    ln -s .singularity.d/actions/test .test && \
    ln -s .singularity.d/runscript singularity

ENV PYTHONIOENCODING=utf-8
WORKDIR ${HOME}

CMD ["bash"]
