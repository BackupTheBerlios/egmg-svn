#!/bin/bash

export CXX='icpc'
export LD='icpc'
export CXXFLAGS='-O3 -ansi -Wall -w2'

gmake -e $@
