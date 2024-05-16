#!/bin/sh
#-*-sh-*-

#
# Copyright © 2009 CNRS
# Copyright © 2009-2020 Inria.  All rights reserved.
# Copyright © 2009 Université Bordeaux
# Copyright © 2014 Cisco Systems, Inc.  All rights reserved.
# See COPYING in top-level directory.
#

HWLOC_top_srcdir="/home/nagp/exageostat-dev/dependinces_installations/hwloc-2.4.0"
HWLOC_top_builddir="/home/nagp/exageostat-dev/dependinces_installations/hwloc-2.4.0"
srcdir="$HWLOC_top_srcdir/utils/hwloc"
builddir="$HWLOC_top_builddir/utils/hwloc"
info="$builddir/hwloc-info"
linuxdir="$HWLOC_top_srcdir/tests/hwloc/linux"

HWLOC_PLUGINS_PATH=${HWLOC_top_builddir}/hwloc/.libs
export HWLOC_PLUGINS_PATH

HWLOC_DEBUG_CHECK=1
export HWLOC_DEBUG_CHECK

HWLOC_DONT_ADD_VERSION_INFO=1
export HWLOC_DONT_ADD_VERSION_INFO

: ${TMPDIR=/tmp}
{
  tmp=`
    (umask 077 && mktemp -d "$TMPDIR/fooXXXXXX") 2>/dev/null
  ` &&
  test -n "$tmp" && test -d "$tmp"
} || {
  tmp=$TMPDIR/foo$$-$RANDOM
  (umask 077 && mkdir "$tmp")
} || exit $?
file="$tmp/test-hwloc-info.output"

set -e
(
  echo "# (default)"
  $info --if synthetic --input "node:2 core:3 pu:4"
  echo
  echo "# --topology"
  $info --if synthetic --input "node:2 core:3 pu:4" --topology
  echo
  echo "# --support"
  $info --if synthetic --input "node:2 core:3 pu:4" --support
  echo
  echo "# --objects"
  $info --if synthetic --input "node:2 core:3 pu:4" --objects
  echo

  echo "# Core range"
  $info --if synthetic --input "node:2 core:3 pu:4" core:2-4
  echo

  echo "# all ancestors of PU range"
  $info --if synthetic --input "node:2 core:3 pu:4" -n --ancestors pu:10-11
  echo
  echo "# Core ancestors of PU range"
  $info --if synthetic --input "node:2 core:3 pu:4" --ancestor core pu:7-9
  echo
  echo "# L2 ancestor of PU"
  $info --if synthetic --input "node:2 core:2 l2:2 l1d:2 pu:2" --ancestor l2 pu:12
  echo
  echo "# L1 ancestor of PU range"
  $info --if synthetic --input "node:2 core:2 l2:2 l1d:2 pu:2" --ancestor l1 -s pu:7-10
  echo

  echo "# Children of L2 and Core of Node, silent"
  $info --if synthetic --input "node:2 core:2 l2:2 l1d:2 pu:2" --children -s l2:1 node:1.core:1
  echo
  echo "# L1d descendants of Core range, silent"
  $info --if synthetic --input "node:2 core:2 l2:2 l1d:2 pu:2" --descendants l1d -s core:1-2
  echo

  echo "# 2 local memory for one PU"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory pu:8
  echo
  echo "# 2 local-or-larger memories for one PU, silent"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory-flags larger -s pu:8
  echo
  echo "# no local-or-larger memory for root, silent"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory-flags larger -s root
  echo
  echo "# no local-or-smaller memory for one PU, silent"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory-flags smaller -s pu:8
  echo
  echo "# 3 local-or-smaller memories for on Package, silent"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory-flags smaller -s pack:1
  echo
  echo "# no strict-local memory for one PU, silent"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory-flags none -s pu:8
  echo
  echo "# 1 strict-local memory for one NUMANode, silent"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory-flags none -s node:1
  echo
  echo "# 12 local-all memories for one PU, silent"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory-flags all\$ -s pu:3
  echo
  echo "# only the smallest locality among 2 local-or-larger memories for one PU, silent"
  $info --if synthetic --input "pack:4 [numa] l3:2 [numa] core:4 pu:2" --local-memory-flags larger --best-memattr locality -s pu:63
  echo
  echo "# only the highest capacity among 2 local-or-larger memories for one PU, silent"
  $info --if synthetic --input "pack:4 [numa(memory=1000000)] l3:2 [numa(memory=1000)] core:4 pu:2" --local-memory-flags larger --best-memattr capacity -s pu:63
  echo

  echo
  echo "# cpukinds for the entire machine"
  $info --if xml --input $linuxdir/fakeheterocpunuma.output all
  echo
  echo "# cpukind for a single PU"
  $info --if xml --input $linuxdir/fakeheterocpunuma.output pu:15
  echo

) > "$file"
/usr/bin/diff -u -w $srcdir/test-hwloc-info.output "$file"
rm -rf "$tmp"
