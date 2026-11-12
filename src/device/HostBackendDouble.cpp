// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.

// The device kernels on the CPU in double precision (see HostBackendImpl.h)
#define IPG_HOST_REAL double
#define IPG_HOST_FACTORY createHostBackendDouble
#include "HostBackendImpl.h"
