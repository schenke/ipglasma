// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_NUCLEUSROLE_H_
#define SRC_NUCLEUSROLE_H_

/**
 * Which of the two colliding nuclei a lattice/nucleon-sampling operation
 * applies to.
 */
enum class NucleusRole {
    /// The first (\f$A\f$) nucleus.
    Projectile,
    /// The second (\f$B\f$) nucleus.
    Target,
};

#endif  // SRC_NUCLEUSROLE_H_
