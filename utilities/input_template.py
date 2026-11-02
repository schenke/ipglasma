# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

"""IP-Glasma input files from a template, as the validation scripts write
them, and the check that a cached run can be reused.

The scripts cache each run in a directory named after the settings they
vary (e.g. the lattice size and the color-charge seed). Other options, such
as the template, the impact parameter or the box size, change the input of
the run without changing that name, so a cached result is only reused if
the input file in its directory gives the same parameters as the input the
run would be written with now (cached_input_mismatch).
"""

import os

# Keys whose values are file or directory paths. The validation scripts write
# them as absolute paths into --ipglasma-path, which can differ between runs
# with the same input (another build directory, --analyze-only from
# elsewhere), so only the name of the file or directory is compared.
PATH_KEYS = ("nucleusQsTableFileName", "nuclearConfigurationsPath",
             "wilsonLinePath")


def input_lines(template, overrides):
    """The lines of the template with the values of the keys in overrides
    replaced; a KeyError for keys that are not in the template."""
    remaining = dict(overrides)
    lines = []
    with open(template) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 2 and tokens[0] in remaining:
                lines.append(f"{tokens[0]} {remaining.pop(tokens[0])}\n")
            else:
                lines.append(line)
    if remaining:
        raise KeyError(f"keys not in the template: {sorted(remaining)}")
    return lines


def write_input(template, overrides, path):
    """Write the template with the values of the keys in overrides
    replaced to path."""
    lines = input_lines(template, overrides)
    with open(path, "w") as f:
        f.writelines(lines)


def parameters(lines):
    """The parameters of input-file lines, key -> value, as IP-Glasma reads
    them: without comments and blank lines, up to EndOfFile."""
    out = {}
    for line in lines:
        tokens = line.split("#", 1)[0].split()
        if tokens == ["EndOfFile"]:
            break
        if tokens:
            out[tokens[0]] = " ".join(tokens[1:])
    return out


def comparable(params):
    """The parameters with the paths of PATH_KEYS reduced to the name of the
    file or directory."""
    return {key: os.path.basename(os.path.normpath(value))
            if key in PATH_KEYS and value else value
            for key, value in params.items()}


def cached_input_mismatch(workdir, template, overrides, name="input",
                          shown=6):
    """None if the input file <workdir>/<name> of a cached run has the
    parameters that the template and the overrides give (comments, the
    order of the lines and the directories of the paths in PATH_KEYS do not
    matter), else the reason why the run cannot be reused, with up to
    `shown` of the parameters that differ."""
    path = os.path.join(workdir, name)
    if not os.path.exists(path):
        return (f"the cached run in {workdir} has no input file to compare "
                "with the current options; remove the directory to rerun it")
    with open(path) as f:
        cached = parameters(f)
    expected = parameters(input_lines(template, overrides))
    same_cached, same_expected = comparable(cached), comparable(expected)
    differ = [f"{key} {cached.get(key, '(missing)')} -> "
              f"{expected.get(key, '(missing)')}"
              for key in sorted(set(cached) | set(expected))
              if same_cached.get(key) != same_expected.get(key)]
    if not differ:
        return None
    more = f", and {len(differ) - shown} more" if len(differ) > shown else ""
    return (f"the cached run in {workdir} was made with other input "
            f"parameters ({', '.join(differ[:shown])}{more}); rerun with "
            "the options of that run, use another --datadir, or remove the "
            "directory to rerun it")
