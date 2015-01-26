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

# script to make a number of setups
# see end of script for example
import sys
import os,errno
import re
from tempfile import mkstemp
from shutil import move,copy
from os import remove, close
import datetime
from numpy import sqrt,array,arange,zeros
import getpass

# forced recursive makedir     
# See http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
# Contributed by user TZOT
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
myusername=getpass.getuser()

# indicate where several files should be located        
codehead='/project/ch4/'+myusername+'/'
codesdir='/freezes'
templatedir='/users/'+myusername+'/hymet_idealized/templates'
headdir='/users/'+myusername+'/csim/fromtempl'
globallog='/users/'+myusername+'/csim/globallog'
eps=1e-8 # A small number (for calculation purposes)
logtext=raw_input('Type something for the log. \n')

# Hardcoded COSMO 1 levels
c1levs=[  22000.00, 21187.52, 20396.25, 19625.89, 18876.13,
          18146.68, 17437.22, 16747.45, 16077.08, 15425.80,
          14793.30, 14179.29, 13583.47, 13005.53, 12445.19,
          11902.13, 11376.06, 10866.69, 10373.71,  9896.84,
           9435.77,  8990.21,  8559.87,  8144.45,  7743.66,
           7357.21,  6984.80,  6626.15,  6280.97,  5948.96,
           5629.83,  5323.31,  5029.09,  4746.91,  4476.46,
           4217.46,  3969.64,  3732.70,  3506.37,  3290.36,
           3084.40,  2888.20,  2701.48,  2523.96,  2355.38,
           2195.45,  2043.90,  1900.45,  1764.84,  1636.78,
           1516.02,  1402.27,  1295.28,  1194.77,  1100.48,
           1012.14,   929.49,   852.28,   780.23,   713.09,
            650.61,   592.52,   538.57,   488.52,   442.11,
            399.09,   359.21,   322.24,   287.93,   256.04,
            226.33,   198.57,   172.54,   148.00,   124.73,
            102.51,    81.13,    60.39,    40.07,    20.00,
              0.00]

# get the length of the large scale forcing file, as a check (since currently no interpolation is used as for e.g. atmos.input)
def getlsflines(lsffile):
    file = open(lsffile,'rb')
    count=0
    timecount=0
    for line in file:
        if len(line.strip())>0:
            if line.strip()[0]!='#':
                count+=1
            elif line.strip()[0]=='#':
                timecount+=1
    return int(count/timecount)

# get maximum height from atmosfile
def getzmax(atmosfile):
    file = open(atmosfile,'rb')
    for line in file:
        if len(line.split())>1:
           lastline = line.split()
    return float(lastline[1])
   
# Root finding: Ridder's algorithms (stable when root is bracketed)
# Used to calculate the distribution of levels
def ridder(f,a,b,tol=1.0e-9):   
    fa = f(a)
    if fa == 0.0: return a
    fb = f(b)
    if fb == 0.0: return b
    if fa*fb > 0.0: print 'Root is not bracketed'
    for i in range(30):
        # Compute the improved root x from Ridder's formula
        c = 0.5*(a + b); fc = f(c)
        s = sqrt(fc**2 - fa*fb)
        if s == 0.0: return None
        dx = (c - a)*fc/s
        if (fa - fb) < 0.0: dx = -dx
        x = c + dx; fx = f(x)
        # Test for convergence
        if i > 0:
            if abs(x - xOld) < tol*max(abs(x),1.0): return x
        xOld = x
        # Re-bracket the root as tightly as possible
        if fc*fx > 0.0:
            if fa*fx < 0.0: b = x; fb = fx
            else:           a = x; fa = fx
        else:
            a = c; b = x; fa = fc; fb = fx
    return None
    print 'Too many iterations'

# This routine does the actual substitutions in the namelist
# for the environmental vars (the ones starting with $) a special flag (lcip) is used, 
# because $ is a special symbol in Python regular expressions
# One of the variables (zlev), is indicated with an lright flag,
# because this variable is used in multiple namelists
# Here the places where substitution are needed
# are indicated on the right hand side of the equality sign
def replaceAll(file_path, searchExp, replaceExp,lright=False,lcip=False):
    #Create temp file
    count=0
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'wb')
    old_file = open(file_path,'rb')
    lremovetrailing=False
    for line in old_file:
        # do not consider trailing lines that look like multi-line variables
        if (lremovetrailing==True) and not(line.strip()=='/' or '=' in line or line.strip()==''):
            continue
        elif (searchExp==line.lstrip()[0:len(searchExp)] or (lright and searchExp in line)):
            # subsitute right hand side
            if(lright):
                sub,plus=re.subn(' *= *'+searchExp+'(.*)($)',' = '+replaceExp+'\r',line)
            # subsitute environmental variable
            elif(lcip):
                sub,plus=re.subn(searchExp[1:]+' *=(.*)($)',searchExp[1:]+' = '+replaceExp+'\r',line)
            # subsitute left hand side
            else:
                sub,plus=re.subn(searchExp+' *=(.*)($)',searchExp+' = '+replaceExp+'\r',line)
            new_file.write(sub)
            count+=plus
            # remove following lines if they do not read solely a / and do not contain =
            if plus>0:
                lremovetrailing=True
        else:
            new_file.write(line)
            lremovetrailing=False
    # print a warning if expression is not replaced a single time
    if(not(lright) and count!=1):
        print "warning: replaced "+searchExp+" by "+replaceExp+" in "+file_path+" "+str(count)+" times" 
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

# This routine adds commas or ; at the end of each line in the cip file
# Needed for e.g. levels            
def convertkeys(dictionary,cip=False):
    if(cip):
        return dict((str(k), str(v)+';') 
            for k, v in dictionary.items())        
    else:
        return dict((str(k), str(v)+',') 
            for k, v in dictionary.items())     

# Current date as string, used in file naming
def todaystr():
    return datetime.date.today().strftime('%Y%m%d')

# Makes exponentially distributed levels
def makelevs(dzmin,zmax,nlevs):
    f = lambda x: zmax-dzmin*((1.0-x**(nlevs-1))/(1.0-x))
    r=ridder(f,1.000000000000001,2.0)
    #print 'growth factor equals '+str(r)
    arr=[0.0]
    for i in range(1,nlevs):
        arr=arr+[dzmin*((1.0-r**i)/(1.0-r))]
    arr[nlevs-1]=zmax-eps
    return arr

# Metadata of a 2d hill                       
class hill:
    def __init__(self,x,y,z,hilltype,hillsideradius_y,hill_width_y=eps,hill_width_x=eps,rotangle=0.0):
        self.x=x
        self.y=y
        self.z=z
        self.hilltype=hilltype
        self.hillsideradius_y=hillsideradius_y
        if hill_width_y==eps:
            self.hill_width_y=hillsideradius_y
        else:
            self.hill_width_y=hill_width_y
        self.hill_width_x=hill_width_x
        self.rotangle=rotangle

# Metadata of a case                                             
class case:
    def __init__(self,lx,ly,hours,lxcrm=None,lycrm=None,dxles=100.0,dzles=40.0,llsf=False,lgsp=True,coriolis=False,nrles=1,nrcrm=1,z0_c=0.035,startlat_tot=15.0,fr_land_c=1.0,zles=makelevs(20,22000,179),rdheight=None,rdtausec=30.0,rasofile_t_is_theta=False,dt_factor=1.0,href_oro=0.0,lsynsat=True):
        self.lx=lx
        self.ly=ly
        if(lxcrm==None):
            self.lxcrm=self.lx
        else:
            self.lxcrm=lxcrm
        if(lycrm==None):
            self.lycrm=self.ly
        else:
            self.lycrm=lycrm            
        self.hours=hours
        self.dxles=dxles
        self.dzles=dzles
        self.rdheight=rdheight # rayleigh damping layer properties 
        self.rdtausec=rdtausec
        self.llsf=llsf # large scale forcings
        self.lgsp=lgsp # microphysics
        self.coriolis=coriolis
        self.rasofile_t_is_theta=rasofile_t_is_theta
        self.nrles=nrles
        self.nrcrm=nrcrm
        self.zles=zles
        self.z0_c=z0_c
        self.startlat_tot=startlat_tot
        self.fr_land_c=fr_land_c
        self.dt_factor=dt_factor
        self.href_oro=href_oro
        self.lsynsat=lsynsat
        self.hills=[]
        self.optlist={} # a dictionary with options
        self.setvals()
    def add2dhill(self,*args,**kwargs):
        self.hills.append(hill(*args,**kwargs))
    def setvals(self):
        self.setflag(self.llsf,'llsf')
        self.setflag(self.lgsp,'lgsp')
        self.setflag(self.coriolis,'lcori')
        self.setflag(self.coriolis,'lcori_deep')
        self.setflag(self.lsynsat,'lsynsat')
        self.setflag(self.rasofile_t_is_theta,'rasofile_t_is_theta')
        self.optlist.update({
        'z0_c':self.z0_c,'fr_land_c':self.fr_land_c,'href_oro':self.href_oro,'hstop':self.hours,'startlat_tot':self.startlat_tot})
    def setflag(self,flag,flagname):
        if flag:
           self.optlist[flagname]='.true.'
        else:
           self.optlist[flagname]='.false.'

# Metadata of a configuration                                                                         
class conf:
    def __init__(self,turb='raschendorfer',res='lores',conv='shall',code='rttov_plus_diag',l2dim=False,dt_factor=1.0,nbounds=3,cloud_num=5.0e08,xkd=0.10):
        self.res=res
        self.l2dim=l2dim
        self.turb=turb
        self.conv=conv
        self.code=code
        self.xkd=xkd
        self.nbounds=nbounds
        self.cloud_num=cloud_num
        self.dt_factor=dt_factor
        self.optlist={}
        self.setvals()
    def setvals(self):
        turbdict={}
        # Dictionary for turbulence options
        turbdict['raschendorfer']={
        'tkhmin':'0.4','lhordiff':'.true.' ,'itype_turb':'3' ,'lprog_tke':'.false.','l3dturb':'.false.','icldm_turb':'2',
        'lkhdef3d':'.false.','lisotrop':'.false.','itype_sher':'0'}
        turbdict['raschendorfer2']={
        'tkhmin':'0.4','lhordiff':'.true.' ,'itype_turb':'3' ,'lprog_tke':'.false.','l3dturb':'.false.','icldm_turb':'2',
        'lkhdef3d':'.false.','lisotrop':'.false.','itype_sher':'2'}
        turbdict['raschendorfer3']={
        'tkhmin':'0.4','lhordiff':'.true.' ,'itype_turb':'3' ,'lprog_tke':'.false.','l3dturb':'.false.','icldm_turb':'2',
        'lkhdef3d':'.false.','lisotrop':'.false.','itype_sher':'3'}
        # Note that the lkhdef3d flag only applies to the hybrid schemes.
        turbdict['smag']={
        'tkhmin':'0.0','lhordiff':'.false.','itype_turb':'11','lprog_tke':'.false.','l3dturb':'.true.' ,'icldm_turb':'1',
        'lkhdef3d':'.false.','lisotrop':'.true.' ,'itype_sher':'0'}
        turbdict['tke']={
        'tkhmin':'0.0','lhordiff':'.false.','itype_turb':'8' ,'lprog_tke':'.true.' ,'l3dturb':'.true.' ,'icldm_turb':'1',
        'lkhdef3d':'.false.','lisotrop':'.true.' ,'itype_sher':'0'}
        turbdict['hybrid']={
        'tkhmin':'0.0','lhordiff':'.false.','itype_turb':'3','lprog_tke':'.false.','l3dturb':'.true.' ,'icldm_turb':'2',
        'lkhdef3d':'.true.' ,'lisotrop':'.true.' ,'itype_sher':'0'}
        turbdict['hybrid2ddef']={
        'tkhmin':'0.0','lhordiff':'.false.','itype_turb':'3','lprog_tke':'.false.','l3dturb':'.true.' ,'icldm_turb':'2',
        'lkhdef3d':'.false.' ,'lisotrop':'.true.' ,'itype_sher':'0'}
        turbdict['tkeadv']={
        'tkhmin':'0.4','lhordiff':'.true.' ,'itype_turb':'3' ,'lprog_tke':'.true.' ,'l3dturb':'.false.','icldm_turb':'2',
        'lkhdef3d':'.false.','lisotrop':'.false.','itype_sher':'0'}
        # Dictionary for convection options
        convdict={}
        convdict['shall']={'nincconv':'10','lconv':'.true.' }
        convdict['noshall']={'nincconv':'1' ,'lconv':'.false.'} 
        convdict['everyt']={'nincconv':'1' ,'lconv':'.true.' }
        # Add dictionaries to configuration
        self.optlist.update(turbdict[self.turb])
        self.optlist.update(convdict[self.conv])
        # Add other options
        self.optlist.update({'cloud_num':self.cloud_num})
        self.optlist.update({'xkd':self.xkd})
        self.setflag(self.l2dim,'l2dim')
        self.setflag(not(self.l2dim),'lperi_y')
    def setflag(self,flag,flagname):
        if flag:
           self.optlist[flagname]='.true.'
        else:
           self.optlist[flagname]='.false.'

# Combine case and configuration to parse into the cip file                                                                                
class confcase:
    def __init__(self,case,conf,casename,confname):
        self.case=case
        self.conf=conf
        self.casename=casename
        self.confname=confname
        self.templatedir=templatedir+'/'+self.casename+'/'
        # a combined dictionary, which includes the dictionaries of case and configuration
        self.optlist={}
        self.optlist.update(self.case.optlist)
        self.optlist.update(self.conf.optlist)       
        self.setdefaults()
        self.parseresolution()
        self.parsehills()
        self.optlist=convertkeys(self.optlist)
        self.optlistcip=convertkeys(self.optlistcip,cip=True)
    def parseresolution(self):
        if self.conf.res=='hires':
            self.dx=self.case.dxles
            self.dy=self.case.dxles        
            self.lx=self.case.lx
            self.ly=self.case.ly
            self.nr=self.case.nrles
            tempz=self.case.zles[::-1]
        elif self.conf.res=='lores':
            self.dx=1000.0
            self.dy=1000.0
            self.lx=self.case.lxcrm
            self.ly=self.case.lycrm
            self.nr=self.case.nrcrm
            tempz=c1levs
        self.vcoordvec=[]
        zmax=getzmax(self.templatedir+'atmos.input')
        for z in tempz:
            if(z<zmax):
                self.vcoordvec.append(z)
        if self.case.rdheight==None:
            self.rdheight=min([0.8*zmax,11500.0,0.8*self.vcoordvec[0]])
        else:
            self.rdheight=self.case.rdheight
        if self.case.rdheight==None:
            self.vcflat=min([0.75*zmax,11357.0,0.75*self.vcoordvec[0]])
        else:
            self.vcflat=self.case.rdheight
        self.dlon=0.00899289*0.001*self.dx
        self.dlat=0.00899289*0.001*self.dy
        almostone=1.0-eps
        self.ie_tot=int(self.lx/self.dx+eps)+2*self.conf.nbounds
        if self.conf.l2dim:
            self.je_tot=1+2*self.conf.nbounds
        else:
            self.je_tot=int(self.ly/self.dy+eps)+2*self.conf.nbounds
        # use quite a few processors (generally good for queueing and efficiency)           
        self.nprocx=int( (self.ie_tot-2*self.conf.nbounds)/16.0+almostone) # make sure each processor gets 16 cells or less
        self.nprocy=int( (self.je_tot-2*self.conf.nbounds)/16.0+almostone) # make sure each processor gets 16 cells or less
        self.dt=self.dx*0.01*self.case.dt_factor*self.conf.dt_factor
        self.nrdtau=int((self.case.rdtausec+eps)/self.dt)
        # vertical levels, get from prelisted???
        self.ke_tot=len(self.vcoordvec)-1
        self.code=codesdir+'/'+self.conf.code
        prefac=1.0/130.0e6 # Efficiency prefactor
        self.walltime=int((prefac*self.case.hours*self.ie_tot*self.je_tot*(self.ke_tot+1)*3600.0)/(self.dt*self.nprocx*self.nprocy)+almostone)
        if self.walltime>24:
           # try again with less cells per core
           prefac=1.0/115.0e6
           self.nprocx=int( (self.ie_tot-2*self.conf.nbounds)/8.0+almostone) # make sure each processor gets 8 cells or less
           self.nprocy=int( (self.je_tot-2*self.conf.nbounds)/8.0+almostone) # make sure each processor gets 8 cells or less
           self.walltime=int((prefac*self.case.hours*self.ie_tot*self.je_tot*(self.ke_tot+1)*3600.0)/(self.dt*self.nprocx*self.nprocy)+almostone)
        if self.walltime>24:
           # try again with less cells per core
           prefac=1.0/100.0e6
           self.nprocx=int( (self.ie_tot-2*self.conf.nbounds)/4.0+almostone) # make sure each processor gets 4 cells or less
           self.nprocy=int( (self.je_tot-2*self.conf.nbounds)/4.0+almostone) # make sure each processor gets 4 cells or less
           self.walltime=int((prefac*self.case.hours*self.ie_tot*self.je_tot*(self.ke_tot+1)*3600.0)/(self.dt*self.nprocx*self.nprocy)+almostone)
        if self.walltime>24:
            sys.exit("Error: the wall time is so large that this job would be utterly inefficient")                                 
        self.resopts={
        'dlon': self.dlon,
        'dlat': self.dlat,
        'ie_tot' : self.ie_tot,
        'je_tot' : self.je_tot,
        'dt' : self.dt,
        'vcoordvec' : ',\r\n'.join(["%16.5f"%i for i in self.vcoordvec]),
        'ke_tot' : self.ke_tot,
        'vcflat' : self.vcflat,
        'rdheight' : self.rdheight,
        'nrdtau' : self.nrdtau,
        }
        self.optlist.update(self.resopts)
        # Generate cross-sections at 8 levels (which ones depends on domain height)
        if(self.vcoordvec[1]>16000):
            self.zlevsparse=array([25.0,125.0,500.0,1000.0,3000.0,6000.0,11000.0,16000.0])+self.case.href_oro
        elif(self.vcoordvec[1]>11000):
            self.zlevsparse=array([25.0,125.0,250.0,500.0 ,1000.0,3000.0,6000.0, 11000.0])+self.case.href_oro
        elif(self.vcoordvec[1]>6000):
            self.zlevsparse=array([25.0,125.0,250.0,500.0 ,1000.0,2000.0, 3000.0, 6000.0])+self.case.href_oro
        elif(self.vcoordvec[1]>3000):
            self.zlevsparse=array([25.0,125.0,250.0,500.0 ,1000.0,1500.0, 2000.0, 3000.0])+self.case.href_oro      
        elif(self.vcoordvec[1]>2000):          
            self.zlevsparse=array([25.0,125.0,250.0,500.0 ,750.0 ,1000.0, 1500.0, 2000.0])+self.case.href_oro         
        elif(self.vcoordvec[1]>1500):          
            self.zlevsparse=array([25.0,125.0,250.0,500.0 ,750.0 ,1000.0, 1250.0, 1500.0])+self.case.href_oro 
        elif(self.vcoordvec[1]>1000):          
            self.zlevsparse=array([25.0,75.0 ,125.0,250.0 ,375.0 ,500.0 , 750.0 , 1000.0])+self.case.href_oro 
        else:
            sys.exit("Error: you have a very shallow domain")
        # special options lists for cip/variables that occur multiple times
        self.optlistright={'zlevhd':',\r\n'.join(["%16.5f"%(i+self.case.href_oro) for i in self.vcoordvec[:1:-1]])+',',
        'zlevsparse':',\r\n'.join(["%16.5f"%i for i in self.zlevsparse])+',',
        'hcomb1':'0.0,'+str(self.case.hours)+',0.1,',
        'hcomb2':'0.0,'+str(self.case.hours)+',0.03333333333333,',
        'hcomb3':'0.0,'+str(self.case.hours)+',0.0416666666666,',
        'hcomb4':'0.0,'+str(self.case.hours)+',0.1,'}
        self.optlistcip={
        '$glocal_code' : '\"'+codehead+'\"',
        '$gcode_vers' : '\"'+self.code+'\"',
        '$gjob_nprocx' : self.nprocx,
        '$gjob_nprocy' : self.nprocy,
        '$gjob_walltime' : '\"%02d:00:00\"'%self.walltime,
        }
    def parsehills(self):
        for hill in self.case.hills:
            self.hillopts['itype_topo']='1'
            self.hillopts['l3dturb_metr']='.true.'         
            self.hillopts['lhill']+=','+'.true.'
            self.hillopts['lhill_2d']+=','+'.true.'
            self.hillopts['hill_type']+=',\''+hill.hilltype+'\''
            self.hillopts['hill_i']+=','+str(hill.x/self.dx+self.conf.nbounds+0.5)
            self.hillopts['hill_j']+=','+str(hill.y/self.dy+self.conf.nbounds+0.5)
            self.hillopts['hillheight']+=','+str(hill.z)
            self.hillopts['hill_rotangle']+=','+str(hill.rotangle)
            self.hillopts['zhillcutfact']+=','+'0.00'
            self.hillopts['hill_combineaction']+=','+'1'
            self.hillopts['hill_width_x']+=','+str(hill.hill_width_x)
            self.hillopts['hill_width_y']+=','+str(hill.hill_width_y)
            self.hillopts['hillsideradius_y']+=','+str(hill.hillsideradius_y)
            self.hillopts['hillasym_x']+=','+'1.0'
            self.hillopts['hillasym_y']+=','+'1.0'
        self.optlist.update(self.hillopts)  
    def setdefaults(self):
        # Include a dummy hill. This makes it easier to parse the hills in a simple way. 
        self.hillopts={
        'itype_topo':'0',
        'l3dturb_metr':'.false.',
        'lhill' : '.false.',
        'lhill_2d' : '.true.',
        'hill_type' : '\'gauss\'',
        'hill_i' : '0.0',
        'hill_j' : '0.0',
        'hillheight' : '0.0',
        'hill_rotangle' : '0.0',
        'zhillcutfact' : '0.00',
        'hill_combineaction' : '1',
        'hill_width_x' : '0.0',
        'hill_width_y' : '0.0',
        'hillsideradius_y' : '0.0',
        'hillasym_x' : '1.0',
        'hillasym_y' : '1.0',
        }
    def copyfiles(self):
        for exp in range(self.nr):
            destdir=headdir+'/'+self.casename+'/'+self.confname+'/'+todaystr()+'exp_%03d'%exp
            print destdir
            mkdir_p(destdir)
            copy(self.templatedir+'atmos.input',destdir)
            copy(self.templatedir+'fluxes.input',destdir)
            copy(self.templatedir+'soil.input',destdir)            
            copy(templatedir+'/cip_globals',destdir)
            if os.path.exists(templatedir+'/rtcoef_msg_1_seviri.dat'):
                copy(templatedir+'/rtcoef_msg_1_seviri.dat',destdir)
                copy(templatedir+'/rtcoef_meteosat_7_mviri.dat',destdir)
            else:
                sys.exit("Error: no RTTOV data present, see README")
            # write to log twice
            writelog(logtext,destdir+'/log')           
            writelog(logtext+'\n destdir='+destdir,globallog)           
            if(self.casename=='arm' and self.case.llsf):
                armforcing(self.vcoordvec,destdir+'/lsf.input')             
            elif((self.casename in ['bomex','bomexhires','bomexextend','bomexdamp'] or 'bomexhismall' in self.casename) and self.case.llsf):
                bomexforcing(self.vcoordvec,destdir+'/lsf.input')
            elif(self.case.llsf):
                sys.exit("Error: no large-scale forcing defined")
    def makenamoptions(self):
        for exp in range(self.nr):
            destdir=headdir+'/'+self.casename+'/'+self.confname+'/'+todaystr()+'exp_%03d'%exp
            for key in self.optlist:
                replaceAll(destdir+'/cip_globals',key,self.optlist[key])
            for key in self.optlistright:
                replaceAll(destdir+'/cip_globals',key,self.optlistright[key],lright=True)
            for key in self.optlistcip:
                replaceAll(destdir+'/cip_globals',key,self.optlistcip[key],lcip=True)
            replaceAll(destdir+'/cip_globals','iseed_noise_t',str(606+37*exp)+',')
            self.dochecks(destdir)
    def dochecks(self,destdir):
        #check consistency: length of lsf.input file with number of vertical levels
        noproblem=True
        if(self.case.llsf):
            noproblem=(self.ke_tot+1==getlsflines(destdir+'/lsf.input'))
        if(noproblem==False):
            print 'length of lsf file does not appear right'
            print 'destdir='+destdir+'/lsf.input'
            print 'ke.tot='+str(self.ke_tot)
            print 'lines in lsf files='+str(getlsflines(destdir+'/lsf.input'))
            sys.exit("Error: problem with large-scale forcings")

# make the large scale forcing for the bomex case
def bomexforcing(zz,outputfile):
    outfile = open(outputfile,'wb')
    lenz=len(zz)    
    tmat=[0,28800]
    lent=len(tmat)
    thf=-2.135e-5
    qtf=-1.2e-8
    ugeo=-8.75
    wsubs=-0.0065
    
    zugeofrac=zeros(lenz)
    zsubsfrac=zeros(lenz)
    zcoolfrac=zeros(lenz)
    zqtadfrac=zeros(lenz)
       
    for j in range(lenz):
        if(zz[j]<700):
            zugeofrac[j]=1
        else:
            zugeofrac[j]=1-1.8e-3*(zz[j]-700.)/8.75
    
    for j in range(lenz):
        if(zz[j]<1500):
            zsubsfrac[j]=zz[j]/1500.
        elif(zz[j]<2100):
            zsubsfrac[j]=1.-(zz[j]-1500.)/(2100.-1500.)
    
    for j in range(lenz):
        if(zz[j]<1500):
            zcoolfrac[j]=1.
        elif(zz[j]<2500):
            zcoolfrac[j]=1.-(zz[j]-1500.)/(2500.-1500.)
    
    for j in range(lenz):
        if(zz[j]<300):
            zqtadfrac[j]=1
        elif(zz[j]<500):
            zqtadfrac[j]=1.-(zz[j]-300.)/(500.-300.)
        
    for i in range(lent):
        outfile.write("# %i \n" %(tmat[i]))
        for j in range(lenz):
            outfile.write("%9.4e %9.4e %5.3f %5.3f %9.4e \n" %(thf*zcoolfrac[j],qtf*zqtadfrac[j],ugeo*zugeofrac[j],0.,wsubs*zsubsfrac[j]))
        if(i<(len(tmat)-1)):
            outfile.write("\n")
    outfile.close()

# make the large scale forcing for the arm case
def armforcing(zz,outputfile):
    outfile = open(outputfile,'wb')
    tmat=[0,3*3600,6*3600,9*3600,12*3600,14.5*3600]
    qtfmat=0.001*array([0.08,0.02,-0.04,-0.1,-0.16,-0.3])/3600.
    thfmat=array([-0.125,0.,0.,-0.08,-0.16,-0.26])/3600.
    lenz=len(zz)
    lent=len(tmat)   
    zfrac=zeros(lenz)
    for j in range(len(zz)):
        if(zz[j]<1000):
            zfrac[j]=1
        elif(zz[j]<2000):
            zfrac[j]=(2000.-zz[j])/1000.
        else:
            zfrac[j]=0.
        
    for i in range(lent):
        outfile.write("# %i \n" %(tmat[i]))
        for j in range(lenz):
            outfile.write("%9.4e %9.4e %5.3f %5.3f %5.3f \n" %(thfmat[i]*zfrac[j],qtfmat[i]*zfrac[j],10.,0.,0.))
        if(i<(lent-1)):
            outfile.write("\n")
    outfile.close()

def writelog(inputtext,outputfile):
    if os.path.exists(outputfile):
        outfile = open(outputfile,'a')
    else:
        outfile = open(outputfile,'wb')
    outfile.write(inputtext)
    outfile.write("\n")
    outfile.write("This case was created "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    outfile.write("\n")
    outfile.write("\n")
    outfile.close()
                                                       
cases={}
cases['bomex']=case(6400.0,6400.0,8.0,lxcrm=64000.0,lycrm=64000.0,dxles=100.0,zles=arange(0.0,3100.0,40.0),llsf=True,lgsp=False,coriolis=True,startlat_tot=15.0,nrles=1,nrcrm=1,z0_c = 0.00012,fr_land_c=0.0,rdheight=2600,lsynsat=False)  
cases['bomexhires']=case(6400.0,6400.0,8.0,lxcrm=64000.0,lycrm=64000.0,dxles=40.0,dt_factor=0.5,zles=arange(0.0,3100.0,40.0),llsf=True,lgsp=False,coriolis=True,startlat_tot=15.0,nrles=1,nrcrm=1,z0_c = 0.00012,fr_land_c=0.0,rdheight=2600)                  
cases['arm']=case(66.67*96,66.67*96,14.0,lxcrm=64000.0,lycrm=64000.0,dxles=66.7,zles=arange(0.0,5500.0,50.0),llsf=True,lgsp=False,coriolis=True,startlat_tot=36.0,nrles=1,nrcrm=1,rdheight=3500)                  
cases['schmidli']=case(40000.0,40000.0,5.0,dxles=100.0,lgsp=False,nrles=1,nrcrm=1,zles=makelevs(10.0,11434.5,179),z0_c = 0.02,rasofile_t_is_theta = True,dt_factor=0.5,rdheight=5500.0)                  
cases['schmidlihires']=case(40000.0,40000.0,5.0,dxles=50.0,lgsp=False,nrles=1,nrcrm=1,zles=makelevs(10.0,11434.5,179),z0_c = 0.02,rasofile_t_is_theta = True,dt_factor=0.5,rdheight=5500.0)                  
cases['kirshbaum']=case(320000.0,60000.0,14.0,dxles=200.0,nrles=1,nrcrm=1,zles=makelevs(10.0,21500.0,179))  
cases['lindaRH85small']=case(64000.0,64000.0,14.0,dxles=200.0,nrles=1,nrcrm=1,zles=makelevs(20.0,21500.0,179),href_oro=500.0)                  
cases['lindaRH85']=case(256000.0,256000.0,14.0,dxles=200.0,zles=makelevs(20,21500.0,179),href_oro=500.0)                  
cases['bomexhighdomain']=case(6400.0,6400.0,8.0,lxcrm=64000.0,lycrm=64000.0,dxles=100.0,zles=makelevs(40.0,10000.0,179),llsf=True,lgsp=False,coriolis=True,nrles=1,nrcrm=1,z0_c = 0.00012,fr_land_c=0.0,rdheight=5500)
cases['grabowski2006']=case(60000.0,60000.0,8.0,dxles=200.0,zles=makelevs(20.0,24000.0,179),coriolis=False,nrles=1,nrcrm=1,z0_c = 0.035,rdheight=16000,dt_factor=0.2)

# Add hills to cases
cases['schmidli'].add2dhill(10000.0,20000.0,1500.0,'cos-plateau',9000.,hill_width_x=500.,hill_width_y=10000000.)
cases['schmidli'].add2dhill(30000.0,20000.0,1500.0,'cos-plateau',9000.,hill_width_x=500.,hill_width_y=10000000.)
cases['schmidlihires'].add2dhill(10000.0,20000.0,1500.0,'cos-plateau',9000.,hill_width_x=500.,hill_width_y=10000000.)
cases['schmidlihires'].add2dhill(30000.0,20000.0,1500.0,'cos-plateau',9000.,hill_width_x=500.,hill_width_y=10000000.)
cases['kirshbaum'].add2dhill(160000.0,30000.0,1500.0,'gauss',20000.,hill_width_x=10000000.,hill_width_y=0.01,rotangle=90.0)

# add the configurations
confs={}
confs['cosmo_smag']=conf(turb='smag',res='hires',conv='noshall')
confs['cosmo_tke']=conf(turb='tke',res='hires',conv='noshall')
confs['c1']=conf()
confs['c1_tkeadvect']=conf(turb='tkeadv')
confs['c1_noshall']=conf(conv='noshall')
confs['c1_hybrid_noshall']=conf(turb='hybrid',conv='noshall')
confs['c1_tke_noshall']=conf(turb='tke',conv='noshall')
confs['c1_smag_noshall']=conf(turb='smag',conv='noshall')
confs['c1_2D']=conf(l2dim=True)
confs['c1_hybrid']=conf(turb='hybrid')
confs['c1_smag']=conf(turb='smag')
confs['simple_shcu']=conf(code='simple_shcu')
confs['simple_shcu_frac']=conf(code='simple_shcu_frac')
confs['cosmo_smag_sandie']=conf(code='mix_sandie',turb='smag',res='hires',conv='noshall')
confs['cosmo_tke_3ddiv']=conf(code='3ddiv',turb='tke',res='hires',conv='noshall',xkd=0.01)

mycases=['bomex']
myconfs=['c1']

# run the whole script
def runme():
    caseconfs={}
    for case in mycases:
        for conf in myconfs:
            #use object composition to make the cases
            caseconfs[case+conf]=confcase(cases[case],confs[conf],case,conf)
            caseconfs[case+conf].copyfiles()          
            caseconfs[case+conf].makenamoptions()

# this statement serves merely the purpose of being able to import 
# classes and routines without running the code            
if __name__=="__main__":
    runme()
