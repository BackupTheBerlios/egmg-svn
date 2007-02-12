#!/bin/sh

export INTEL_VERSION=90
export INTEL_PATH=/opt/intel/compiler90/

bjam -sTOOLS=intel-linux $*
