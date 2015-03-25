**LICENSE**

*These programs/scripts are free software: you can redistribute it and/or modify*
*it under the terms of the GNU Lesser General Public License as published by*
*the Free Software Foundation, either version 3 of the License, or*
*(at your option) any later version.*

*This program is distributed in the hope that it will be useful,*
*but WITHOUT ANY WARRANTY; without even the implied warranty of*
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*
*GNU General Public License for more details.*

*You should have received a copy of the GNU Lesser General Public License*
*along with this program.  If not, see <http://www.gnu.org/licenses/>.*
*Copyright (C) 2013-2014 Steven Boeing, ETHZ*

**CONTACT**

steven."lastname with oe" (at) env.ethz.ch
OR
sjboing (at) "the usual g-mail suffix"

**IMPORTANT**

* The RTTOV data (rtcoef_msg_1_seviri.dat, rtcoef_meteosat_7_mviri.dat) that are needed to use this package are currently NOT provided with the package! Download these and put them in the templates directory.

**DESCRIPTION**

This is a collection of scripts for setting up and running idealized (including LEM) simulations with COSMO. These scripts by default assume that :
* The hymet_idealized directory is directly under home (directories as at cscs: /users/...)
* Cip is installed (see notes in the cip directory. Some paths may need to be adjusted
  Provided with this package, thanks Oli Fuhrer and Daniel Leuenberger for making this available)
* Simulations are run in a csim/fromtempl directory in the home folder

This directory includes the following scripts
* makesims_user.py: For setting up cases with a certain COSMO configuration
* runsims_user.py: for starting these cases (can also be done by hand using cip)
* neatpostproc_user.py: Generic postprocessing for idealized cases
* variablelist.py: Lists of "dir" (direct output) and "der" (derived) variables to analyze
* areaspectra.py: code for calculating spectra (thanks Juerg Schmidli)
* kirshbaummaker.py: code for setting up a case in the spirit of
  Daniel J. Kirshbaum, 2011: Cloud-Resolving Simulations of Deep Convection over a Heated Mountain. J. Atmos. Sci., 68, 361–378. 
* cosmo_utils.py: code used by kirshbaummaker.py (thanks Juerg Schmidli)

And a templates directory with basic information on the cases.

A modified version of these scripts has been used by Davide Panosetti during his MSc thesis
A manuscript for MetZ is in preparation

** OUTPUT FILES **

* xzlv.*.nc: mean values in x,z plane, variables outputted on model levels
* xzz.*.nc: mean values in x,z plane, variables outputted on height levels
* interp.*.nc: mean values in x,z plane, variables outputted on model levels, interpolated to a number of height levels
* crossxy.variablename.*.nc: cross sections in x,y plane at different levels
* crosslv.*.nc: cross sections in x,z plane, variables outputted on model levels
* crossz.*.nc: cross sections in x,z plane, variables outputted on height levels
* crossinterp.*.nc: cross sections in x,z plane, variables outputted on model levels, interpolated to a number of height levels
* prof1d.*.nc: mean profiles, variables outputted on height levels, interpolated to height levels
* interp1d.*.nc: mean profiles, variables outputted on model levels, interpolated to a number of height levels, values below topography not taken into account
* hovlv.*.nc: hovmoeller diagrams, derived from model level output
* hovz.*.nc: hovmoeller diagrams, derived from height level output
* hovlv.*.nc: domain averages (many variables similar to hovmoeller diagrams), derived from model level output
* hovz.*.nc: domain averages (many variables similar to hovmoeller diagrams), derived from height level output
* clouds.*.tar: tar file with 3d cloud fields (no gz, as cloud fields themselves are already compressed)

**KNOWN ISSUES/FEATURES**

* The postprocessing is currently only written for cases with a ridge-like
  topography (i.e. in the averaging it is assumed one dimension is homogeneous in height).
* In COSMO, what is outputted as TKE currently depends on both the choice of turbulence scheme and
  whether you dump height interpolated or model level output
* For runsims, the correct partition names depend on the system on which you
  do postprocessing (use sinfo to find out more). For very large cases, one could also
  run out of memory (use node with more memory or rewrite script) or processing time
  (use the long queue on e.g. pilatus)
* Installing netcdf4python on daint : we got it running with the following steps (26 Jan 2015):


  module swap PrgEnv-pgi PrgEnv-gnu


  module load zlib


  module load hdf5-parallel/1.8.13_ftn


  module load netcdf-hdf5parallel/4.3.2_ftn


  export NETCDF4_DIR=/apps/pilatus/netcdf-hdf5parallel/4.3.2_ftn/gnu_482


  export HDF5_DIR=/apps/pilatus/hdf5-parallel/1.8.13_ftn/gnu_482


  export CC=mpicc


  export C_INCLUDE_PATH="/apps/pilatus/mvapich2/1.9/gcc-4.8.2/include/"


  python setup.py


  export PYTHONPATH="$PYTHONPATH:/users/pdavide/lib/python2.7/site-packages/"


  python setup.py install --prefix=/users/YOURUSERNAME

  Possibly, you will also need to:


  module load netcdftime


  module load mvapich2/1.9
