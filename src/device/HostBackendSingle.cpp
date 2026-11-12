// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

// The device kernels on the CPU in single precision (see HostBackendImpl.h):
// the arithmetic of a GPU without double precision, e.g. of an Apple M
// processor
#define IPG_HOST_REAL float
#define IPG_HOST_FACTORY createHostBackendSingle
#include "HostBackendImpl.h"
