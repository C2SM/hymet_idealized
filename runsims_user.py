# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# Copyright (C) 2013-2014 Steven Boeing, ETHZ

# script for running a set of cases
# requires cip!
import glob
import os
import getpass
      
myusername=getpass.getuser()

headdir='/users/'+myusername+'/csim/fromtempl'

mycases=['bomex']
myconfs=['c1']
expglob='20150112exp_000' # which experiments to select

def intersect(a, b):
     return list(set(a) & set(b))

for case in mycases:
    for conf in myconfs:
        # find underlying experiments
        curdir=headdir+'/'+case+'/'+conf+'/'
        exps=glob.glob(curdir+expglob)
        subdirs=[curdir+ i for i in os.walk(curdir).next()[1]]
        # make sure experiment corresponds to actual working case
        for exper in intersect(exps,subdirs):
            os.chdir(exper)
            os.system('cip clean')
            os.system('cip start')          
