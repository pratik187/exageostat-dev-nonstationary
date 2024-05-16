#!/bin/bash -x
# StarPU --- Runtime system for heterogeneous multicore architectures.
#
# Copyright (C) 2017-2021  Université de Bordeaux, CNRS (LaBRI UMR 5800), Inria
#
# StarPU is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 2.1 of the License, or (at
# your option) any later version.
#
# StarPU is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Lesser General Public License in COPYING.LGPL for more details.
#
# Test a model-checking program with simgrid model checker

SIMGRID_MC=
abs_top_srcdir=/home/nagp/exageostat-dev/dependinces_installations/starpu-1.3.9
abs_builddir=/home/nagp/exageostat-dev/dependinces_installations/starpu-1.3.9/tests/model-checking

set -e

[ -x "$SIMGRID_MC" ] || exit 77

MC_FLAGS=--cfg=model-check/reduction:none

# makes it much longer actually
#MC_FLAGS+=--cfg=contexts/factory:ucontext
#MC_FLAGS+=--cfg=model-check/sparse-checkpoint:yes
#MC_FLAGS+=--cfg=model-check/visited:1000

test() {
	$SIMGRID_MC $abs_builddir/$1 $abs_top_srcdir/tests/model-checking/platform.xml MAIN $MC_FLAGS
}
