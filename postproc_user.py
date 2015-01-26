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

import glob
import os
import getpass

myusername=getpass.getuser()

headdir='/home/'+myusername+'/csim/fromtempl'
mycases=['schmidlihires']
myconfs=['cosmo_tke_sandie']
exps=['20150112exp_000']

# make a job from a template
def writejob(case,conf,exper,hours,outputfile,partition):
    outfile = open(outputfile,'wb')
    outfile.write("""#!/bin/bash
#SBATCH --job-name="python job" 
#SBATCH --nodes=1
##SBATCH --dependency=afterok:000000 #you can make this job dependent
""")
    outfile.write("#SBATCH --time=%02d:00:00 \n"%hours)
    outfile.write("#SBATCH --partition="+partition+" \n")
    outfile.write("#SBATCH --output=/users/"+myusername+" \n")
    outfile.write("#SBATCH --error=/users/"+myusername+"/hymet_idealized/slurm-ERR.log")
    outfile.write("""#SBATCH --mem=32g
#SBATCH --exclusive

#SBATCH --ntasks-per-node=1
. /etc/profile.d/modules.bash
ulimit -a
""")
    outfile.write("cd /users/"+myusername+"/hymet_idealized\n")
    outfile.write("python neatpostproc_user.py "+case+" "+conf+" "+exper+"\n")
    outfile.close()
    outfile.close()

# define cases by size
bigcases=['kirshbaum','lindaRH85']
middlecases=['lindaRH85small','schmidli','bomexhires','grabowski2006','schmidlihires']
bigconfs=['cosmo_smag','cosmo_tke','cosmo_tke_sandie']

# make temporary jobs for postprocessing, and submit them
i=0    
for case in mycases:
    for conf in myconfs:
        for exper in exps:
            partition='normal'
            if (conf in bigconfs and case in bigcases):
                hours=47
                partition='long'
            elif (conf in bigconfs and case in middlecases):
                hours=23
            elif(case in bigcases or conf in bigconfs):
                hours=6
            else:
                hours=1
            if os.path.isdir('/scratch/daint/'+myusername+'/exe/csim/fromtempl/'+case+'/'+conf+'/'+exper+'/output/'):
                writejob(case,conf,exper,hours,"tempjob.%03d"%i,partition)
                os.system('sbatch tempjob.%03d'%i)
                i=i+1
  
