// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

#ifndef SRC_LATTICEINDEX_H_
#define SRC_LATTICEINDEX_H_

// The memory layout shared by every per-site lattice field and by the FFT
// arrays: site (ix, iy) is stored at ix * N + iy (ix outer, iy inner).

/**
 * Site index of the transverse position (\p ix, \p iy) on an
 * \p N x \p N lattice. This is the one place the memory layout of every
 * per-site field is defined (\p ix outer, \p iy inner); use it, or
 * Lattice::positionFromXY(), instead of writing the formula by hand.
 * \param[in] ix Site index in \f$x\f$, `0 <= ix < N`.
 * \param[in] iy Site index in \f$y\f$, `0 <= iy < N`.
 * \param[in] N Lattice side length.
 * \return The index into the per-site vectors.
 */
inline int latticeIndex(int ix, int iy, int N) { return ix * N + iy; }
/**
 * The \f$x\f$ index of site \p pos; inverse of latticeIndex().
 * \param[in] pos Index into the per-site vectors.
 * \param[in] N Lattice side length.
 * \return The site's \f$x\f$ index.
 */
inline int latticeX(int pos, int N) { return pos / N; }
/**
 * The \f$y\f$ index of site \p pos; inverse of latticeIndex().
 * \param[in] pos Index into the per-site vectors.
 * \param[in] N Lattice side length.
 * \return The site's \f$y\f$ index.
 */
inline int latticeY(int pos, int N) { return pos % N; }

#endif  // SRC_LATTICEINDEX_H_
