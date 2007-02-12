#!/bin/bash

export CXX='CC'
export LD='CC'
export CXXFLAGS='+w2'

gmake -e $@
