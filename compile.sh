#!/bin/bash
set -e -u

module unload cray-netcdf cray-hdf5 fre
module unload PrgEnv-pgi PrgEnv-intel PrgEnv-gnu PrgEnv-cray
module load PrgEnv-intel/6.0.3
module swap intel intel/16.0.3.210

# DO NOT LOAD fre/bronx-16: it changes perl version to 5.28.0, which seem to create problems
# for fpx3. Default /usr/bin/perl (v.5.18.2) apparently works fine.
#module load fre/bronx-16
#module load cray-hdf5/1.8.16

make
