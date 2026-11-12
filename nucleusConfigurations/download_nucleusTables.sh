#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
# This file is part of IP-Glasma; the license text is in LICENSE.

# Downloads the nucleon configurations of nuclear-structure models
# (nucleonPositionsFromFile 1, polarized deuterons) into this folder from
# https://bitbucket.org/wayne_state_nuclear_theory/glaubernucleusconfigs,
# whose README gives their sources. The commit is fixed, and the files are
# checked against SHA256SUMS, so every download gives the same
# configurations. To move to another commit, change COMMIT and regenerate
# SHA256SUMS (sha256sum *.bin.in > SHA256SUMS in a fresh download).
#
# Usage: cd nucleusConfigurations && bash download_nucleusTables.sh
set -euo pipefail
cd "$(dirname "$0")"

REPO_URL="https://bitbucket.org/wayne_state_nuclear_theory/glaubernucleusconfigs.git"
COMMIT="945fe00aabb5cd83af350a799f9646dc138d956d"  # 2026-09-11
TMP="glaubernucleusconfigs"

rm -rf "$TMP"
git init -q "$TMP"
git -C "$TMP" fetch -q --depth 1 "$REPO_URL" "$COMMIT"
git -C "$TMP" checkout -q FETCH_HEAD
mv "$TMP"/nucleusConfigurations/*.bin.in .
mv "$TMP"/polarizedConfigurations/*.bin.in .
rm -rf "$TMP"

if command -v sha256sum > /dev/null; then
    sha256sum -c --quiet SHA256SUMS
else
    shasum -a 256 -c --quiet SHA256SUMS
fi
echo "downloaded and verified the configurations of commit ${COMMIT:0:7}"
