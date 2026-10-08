#! /usr/bin/env python3
"""
     This script fetches an individual IP-Glasma event from the hdf5 database
     written with writeOutputsToHDF5 1 (see OUTPUT.md). It writes the hydro
     file (the input for the MUSIC fluid dynamic simulation) or the text
     T^{mu nu} file of the event at a given proper time.
"""

from sys import argv, exit
import numpy as np
import h5py

PREFIXES = {0: "epsilon-u-Hydro-t", 1: "Tmunu-t"}


def print_help():
    print("{0} database_filename event_id output_type tau".format(argv[0]))
    print("    output_type: 0 for the hydro file, 1 for the text Tmunu file")
    print("    tau: the proper time as written in the file names, e.g. 0.4;")
    print("         without it, the available times are listed")


def header_values(header):
    """Returns the `key= value` pairs of a file header as a dict."""
    tokens = header.replace('#', '').split()
    return {tokens[i][:-1]: tokens[i + 1] for i in range(len(tokens) - 1)
            if tokens[i].endswith('=')}


def available_times(event_group, prefix, event_idx):
    """Returns the proper times of the files of the event with prefix."""
    suffix = "-{0:d}.dat".format(event_idx)
    return sorted((name[len(prefix):-len(suffix)] for name in event_group
                   if name.startswith(prefix) and name.endswith(suffix)),
                  key=float)


def get_dataset(database_path, output_type, time_stamp, event_idx):
    """Returns the dataset and its file name, or exits with an error."""
    hf = h5py.File(database_path, "r")
    event_name = "event-{0:d}".format(event_idx)
    event_group = hf.get(event_name)
    if event_group is None:
        print("{0} has no {1}".format(database_path, event_name))
        exit(1)
    prefix = PREFIXES[output_type]
    file_name = "{0}{1}-{2:d}.dat".format(prefix, time_stamp, event_idx)
    if time_stamp is None or file_name not in event_group:
        times = available_times(event_group, prefix, event_idx)
        print("available times of {0}*-{1:d}.dat: {2}".format(
            prefix, event_idx, " ".join(times) if times else "none"))
        exit(1 if time_stamp is not None else 0)
    return event_group[file_name], file_name


def fetch_an_IPGlasma_event_Tmunu(database_path, time_stamp, event_idx):
    print("fetching an IP-Glasma event Tmunu with event id: {} at tau = {} "
          "fm from {}".format(event_idx, time_stamp, database_path))
    temp_data, file_name = get_dataset(database_path, 1, time_stamp,
                                       event_idx)
    data_header = temp_data.attrs["header"].decode('UTF-8').replace('#', '')
    nx = temp_data.attrs["nx"]
    ny = temp_data.attrs["ny"]
    if len(temp_data) != nx*ny:
        print("{0} has {1} rows, expected {2}".format(
            file_name, len(temp_data), nx*ny))
        exit(1)

    # one line per grid point, iy outer and ix inner
    output_data = np.zeros([len(temp_data), 12])
    output_data[:, 2:] = temp_data
    idx = 0
    for iy in range(ny):
        for ix in range(nx):
            output_data[idx, 0] = ix
            output_data[idx, 1] = iy
            idx += 1
    np.savetxt(file_name, output_data, fmt=('%i  %i' + '  %.6e'*10),
               header=data_header)
    return file_name


def fetch_an_IPGlasma_event(database_path, time_stamp, event_idx):
    print("fetching an IP-Glasma event with event id: {} at tau = {} fm "
          "from {}".format(event_idx, time_stamp, database_path))
    temp_data, file_name = get_dataset(database_path, 0, time_stamp,
                                       event_idx)
    data_header = temp_data.attrs["header"].decode('UTF-8').replace('#', '')
    x_size = temp_data.attrs["x_size"]
    y_size = temp_data.attrs["y_size"]
    nx = temp_data.attrs["nx"]
    ny = temp_data.attrs["ny"]
    # the grid spacing from the grid size, since the dx and dy in the
    # header are rounded to 6 digits
    dx = x_size/nx
    dy = y_size/ny
    values = header_values(data_header)
    neta = int(values["etamax"])
    deta = float(values["deta"])
    if len(temp_data) != neta*nx*ny:
        print("{0} has {1} rows, expected {2}".format(
            file_name, len(temp_data), neta*nx*ny))
        exit(1)

    # one line per grid point, eta outer, then x, then y; the eta, x and y
    # columns were not stored
    output_data = np.zeros([len(temp_data), 18])
    output_data[:, 3:] = temp_data
    idx = 0
    for ieta in range(neta):
        eta_local = -(neta - 1)/2.*deta + ieta*deta
        for ix in range(nx):
            x_local = -x_size/2. + ix*dx
            for iy in range(ny):
                y_local = -y_size/2. + iy*dy
                output_data[idx, 0] = eta_local
                output_data[idx, 1] = x_local
                output_data[idx, 2] = y_local
                idx += 1
    np.savetxt(file_name, output_data, fmt='  '.join(['%.6e']*18),
               header=data_header)
    return file_name


if __name__ == "__main__":
    try:
        database_filename = str(argv[1])
        event_id = int(argv[2])
        type_flag = int(argv[3])
    except (IndexError, ValueError):
        print_help()
        exit(1)
    if type_flag not in PREFIXES:
        print_help()
        exit(1)
    time_stamp_str = argv[4] if len(argv) > 4 else None

    if type_flag == 0:
        fetch_an_IPGlasma_event(database_filename, time_stamp_str, event_id)
    else:
        fetch_an_IPGlasma_event_Tmunu(database_filename, time_stamp_str,
                                      event_id)
