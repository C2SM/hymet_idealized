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




# A PYTHON SCRIPT TO POSTPROCESS ALL COSMO OUTPUT IN A DIRECTORY AT ONCE
# SEPARATE SCRIPTS SHOULD DO THE PLOTTING OF OUTPUT
# AIMED AT A COMBINATION OF SPEED; READABILITY AND A SOMEWHAT LIMIT MEMORY USAGE
# I.E. ABILITY TO POSTPROCESS LARGE DATA ON A SINGLE NODE

# detect all cosmo output files, their type, and find the corresponding time
# incrementally add the output to an output file
# includes masked/sampled statistics (like traditional LES models)
# masks can also be used to average over pre-selected subdomain (e.g. region around hill)
# but currently such domains are not implemented
# outputs to a file which is interpolated to height levels and saved with f4 precision (for file size limitation)

from numpy import *
import numpy as np

from netCDF4 import Dataset
import os,errno
from time import clock
from scipy import weave #Embed code in C
from scipy.weave import converters
import sys # system functions
import glob # a libary for regular expressions in file names
import shutil # for copying the files to project
import areaspectra as asp # A separate library for computing spectra
import tarfile # For compressing the cloud data
from optparse import OptionParser # Flags to add later: store interpolated output?
from variablelist import * # A separate file contains the actual variable list
import getpass     
myusername=getpass.getuser()

nboundlines=3
     
interpout=True # also store interpolated output on regular grid for fast data exploration
loadmpl=False # Load matplotlib
lbud=True # Try to calculate budget variables?
lzlib=True # Compress output using zlib?
lcross=True # Produce cross-sections? Mais oui
lclouds=True # Produce tar-ball with 3d cloud fields for storage
lsats=True # Produce file with satellite stuff
np.seterr(invalid='ignore') # don't wine about nans

# cosmo constants for calculations
pref=1.0e5 # ref pressure
rd=287.05 # gas constant, dry air
rv=461.51 # gas constant, water vapor
cpd=1005.0 # heat capacity at constant pressure
rlv=2.501e6 # latent heat of condensation
riv=2.835e6 # latent heat of freezing
grav=9.81 # gravitational acceleration

# constants used in calculation of saturation pressure (COSMO specific)
b1=610.78 
b2w=17.2693882
b2i=21.8745584
b3=273.16
b4w=35.86
b4i=7.66

start=clock()

# currently not used, but useful when analyzing problems
if loadmpl:
    import matplotlib
    matplotlib.use('agg') # first define agg output, then continue to load rest of matplotlib and pyplot 
    from matplotlib import *
    from matplotlib.pyplot import *
    
# forced makedir
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# specify different cosmo outputfile types
def make_filelist():
    global filelist
    types = {'level':'lfff*.nc','cloud':'lf*.nca','crossxy':'lf*.ncb','height':'lfff*.ncc','sats':'lfff*.ncs'} # the file types
    filelist={}
    for i in types:
        filelist[i]=glob.glob(fulldir+types[i])
        filelist[i].sort() 
        
# create a tarfile (for cloud data)
def make_tarfile(output_filename, infiles):
    with tarfile.open(output_filename,'w') as tar:
        for infile in infiles:
            tar.add(infile)
                          
# class for masks/conditional sampling
class mask:
    def __init__(self):
        pass
    def setfield(self,field):
        self.field=field

# mask for spectra, currently the same as regular mask class
class specmask(mask):
    pass

# set values in an array of heigths to at least the surface value
# useful for integration of e.g. water vapor anomaly
# (whenever dz needed, cosmo simply interpolates on)
# uses C-code (weave) for speed
def replace_h_hlower(h,hlower,tmax,kmax,jmax,imax):
    code = """
    int i, j, k, t;
    for (t=0; t<tmax; t++) {
        for (k=0; k<kmax; k++) {
           for (j=0; j<jmax; j++) {
              for (i=0; i<imax; i++) {
                 if(h(t,k,j,i)<hlower(j,i)) {
                     h(t,k,j,i)=hlower(j,i);
                 }
              }     
           }
        }
    }
    """
    weave.inline(code,['h','hlower','tmax','kmax','jmax','imax'],type_converters=converters.blitz,compiler='gcc')
 
# weave interpolation to height levels
# takes advantage of known cosmo 2d layout
# uses C-code (weave) for speed
def int_to_height(outarr,outheights,inarr,inheights,tmax,kmax,jmax,kmaxout):
    code = """
    int t,k,j,ks;
    int kstore[tmax][jmax];
    for (t=0; t<tmax; t++) {
        for (j=0; j<jmax; j++) {
            kstore[t][j]=kmaxout-1;
        }
    }
    for (t=0; t<tmax; t++) {
        for (k=0; k<kmax-1; k++) {
           for (j=0; j<jmax; j++) {
                ks=kstore[t][j];
                while(outheights(ks)>inheights(t,k+1,j) and ks>0) {
                    outarr(t,ks,j)=inarr(t,k,j)+(outheights(ks)-inheights(t,k,j))*(inarr(t,k+1,j)-inarr(t,k,j))/(inheights(t,k+1,j)-inheights(t,k,j));
                    ks=ks-1;
                }
                if(ks==0 and outheights(ks)>(inheights(t,k+1,j)-1e-8)) {
                    outarr(t,ks,j)=inarr(t,k,j)+(outheights(ks)-inheights(t,k,j))*(inarr(t,k+1,j)-inarr(t,k,j))/(inheights(t,k+1,j)-inheights(t,k,j));
                }                
                kstore[t][j]=ks;
            }
        }             
    }
    """
    weave.inline(code,['outarr','outheights','inarr','inheights','tmax','kmax','jmax','kmaxout'],type_converters=converters.blitz,compiler='gcc')

                      
# prepare spectral data into 2d arrays, containing stretches
# along all points that match sampling criteria
# currently assumes xz-type geometry
# uses C-code (weave) for speed
def prepare_spectra_xz(dataout,data,mask,tmax,kmax,jmax,imax,nsegmentsmax):
    code = """
    int i, j, k, t, nsegments;
    nsegments=0;
    for (t=0; t<tmax; t++) {
        for (k=0; k<kmax; k++) {
            for (i=0; i<imax; i++) {
                if(mask(t,k,i)) {
                     for (j=0; j<jmax; j++) {
                        dataout(nsegments,j)=data(t,k,j,i);
                     }
                     nsegments=nsegments+1;
                }
            }     
        }
    }
    """
    weave.inline(code,['dataout','data','mask','tmax','kmax','jmax','imax','nsegmentsmax'],type_converters=converters.blitz,compiler='gcc')
         
# mean value across 2d slab
# note: to calculate a true mean for a run with topography, including filled values, you need to correct for the topography
def mean2d(inputfield):
    meanfield=mean(mean(inputfield,axis=-2, dtype=np.float64),axis=-1, dtype=np.float64)
    return meanfield
    
# mean value across 1d line      
def mean_1d(inputfield,geometry):
    if(geometry=='xz'):
       meanfield=mean(inputfield,axis=-2, dtype=np.float64)
    elif(geometry=='yz'):
       meanfield=mean(inputfield,axis=-1, dtype=np.float64)
    return meanfield

# extraction used for cross-sections
def extract_1d(inputfield,geometry):
    if(geometry=='xz'):
       extrfield=inputfield[:,:,0,:]
    elif(geometry=='yz'):
       extrfield=inputfield[:,:,:,0]
    return extrfield
        
# deviations with respect to a 1d line      
def deviation_1d(inputfield,geometry):
    if(geometry=='xz'):
       meanfield=mean(inputfield,axis=2, dtype=np.float64)
       delta=inputfield-meanfield[:,:,None,:]
    elif(geometry=='yz'):
       meanfield=mean(inputfield,axis=3, dtype=np.float64)
       delta=inputfield-meanfield[:,:,:,None]
    return delta

# command to get a single variable from a file
def var_from_file(dataset,key):
    try:
        return dataset.variables[(key)][:,:,nboundlines:-nboundlines,nboundlines:-nboundlines]
    except:
        return(nan)  

# a class to include some methods needed for both helper objects and netcdf objects
class get_variable_class():
    # gv being "get variable"
    def gv(self,key):
        if key in self.varkeys:
            return self.data.variables[(key)][:,:,nboundlines:-nboundlines,nboundlines:-nboundlines]
        else:
            return(self.data.variables[('P')][:,:,nboundlines:-nboundlines,nboundlines:-nboundlines]*nan)
    # gdim being "get dimension"
    def gdim(self,key):
        try:
            return self.data.variables[(key)][:]
        except:
            try:
                return xrange(len(self.data.dimensions[(key)]))
            except:
                return([nan]) 
                
# a class to store derived variables from output which are needed relatively often       
class nchelper(object,get_variable_class):
    def __init__(self,geom):
        self.data=[]
        self.geom=geom
        self.tstep=0         
    def update(self):
        self.tstep=self.tstep+1         
    def calcrhoqc(self):
        p=self.gv('P')
        t=self.gv('T')
        hydrotot=self.gv('QC')+self.gv('QR')+self.gv('QS')+self.gv('QG')+self.gv('QI')
        # if no ice species defined, correct precipitation
        if(isnan(hydrotot[0,0,0,0])):
            hydrotot=self.gv('QC')+self.gv('QR')+self.gv('QS')+self.gv('QI')
        if(isnan(hydrotot[0,0,0,0])):
            hydrotot=self.gv('QC')  
        trho=t*(1+(rv/rd-1)*self.gv('QV')-hydrotot)
        self.rhoh=p/(rd*trho)       
        qci=self.gv('QC')+self.gv('QI')
        # if ice species are not defined, correct qci
        if(isnan(hydrotot[0,0,0,0])):
           qci=seld.gv('QC')
        self.qt=self.gv('QC')+self.gv('QV')+self.gv('QI')
        if(isnan(self.qt[0,0,0,0])):
           self.qt=self.gv('QC')+self.gv('QV')
        iexnf=(p/pref)**(-rd/cpd) #inverse exner function (not stored)
        # buoyancy potenital temperature
        self.thetarhoh=trho*iexnf
        # cloud mask
        self.cld=(qci>1.0e-6)

# derived variables specific to level output             
class nclevelhelper(nchelper):
    def update(self,data):
        self.data=data
        self.varkeys=self.data.variables.keys()  
        self.calcrhoqc()
        h=self.gv('HHL')
        w=self.gv('W')
        dh=h[:,1:,:,:]-h[:,:-1,:,:]
        # integrated paths of water species
        self.vwp=-sum(self.rhoh*self.gv('QV')*dh,axis=1,dtype=float64)
        self.cwp=-sum(self.rhoh*self.gv('QC')*dh,axis=1,dtype=float64)
        self.rwp=-sum(self.rhoh*self.gv('QR')*dh,axis=1,dtype=float64)
        self.gwp=-sum(self.rhoh*self.gv('QG')*dh,axis=1,dtype=float64)
        self.swp=-sum(self.rhoh*self.gv('QS')*dh,axis=1,dtype=float64)
        self.iwp=-sum(self.rhoh*self.gv('QI')*dh,axis=1,dtype=float64)
        self.wmax=nanmax(w,axis=1)
        self.wmin=nanmin(w,axis=1)
        self.sshf=self.gv('HFLUX')[:,-1,:,:]
        self.slhf=self.gv('EFLUX')[:,-1,:,:]
        self.h=h
        self.hh=0.5*(h[:,1:]+h[:,:-1])
        self.wh=0.5*(w[:,1:]+w[:,:-1])
        # cloud and cloud updraft fractions
        self.cldqcifrac=1.0*(sum(self.cld,axis=1)>0)
        self.cldqciw1frac=1.0*(sum(self.cld*(self.wh>1.0),axis=1)>0)
        # wind at lowest level
        self.llwind=sqrt(self.gv('U')[:,-1,:,:]**2+self.gv('V')[:,-1,:,:]**2)
        super(nclevelhelper,self).update()

# derived variables specific to height output                                            
class ncheighthelper(nchelper):
    def update(self,data):
        self.data=data
        self.varkeys=self.data.variables.keys()     
        self.calcrhoqc()
        if(self.tstep==0):
            self.make_topo()
        ### calculate buoyancies, and fix output
        if(self.geom=='xz'):
            mthetarhoh=mean2d(self.thetarhoh*self.topomask[None,:,None,:])[:,:,None,None]/mean2d(self.topomask[None,:,None,:])[:,:,None,None]
            self.buoy=grav*(self.thetarhoh*self.topomask[None,:,None,:]-mthetarhoh)/(mthetarhoh*self.topomask[None,:,None,:])
        elif(self.geom=='yz'):
            mthetarhoh=mean2d(self.thetarhoh*self.topomask[None,:,:,None])[:,:,None,None]/mean2d(self.topomask[None,:,:,None])[:,:,None,None]
            self.buoy=grav*(self.thetarhoh*self.topomask[None,:,:,None]-mthetarhoh)/(mthetarhoh*self.topomask[None,:,:,None])
        whereinf=isinf(self.buoy);
        self.buoy[whereinf] = nan;
        nanbuoy=self.buoy[:,:,:,:]
        wherefin = isfinite(nanbuoy);
        nanbuoy[wherefin==False] = 0.0;
        ### calculate water vapor anomalies, and fix output
        qv=self.gv('QV')
        if(self.geom=='xz'):
            mqv=mean2d(qv*self.topomask[None,:,None,:])[:,:,None,None]/mean2d(self.topomask[None,:,None,:])[:,:,None,None]
            self.qvp=(qv*self.topomask[None,:,None,:]-mqv)        
        elif(self.geom=='yz'):
            mqv=mean2d(qv*self.topomask[None,:,:,None])[:,:,None,None]/mean2d(self.topomask[None,:,:,None])[:,:,None,None]
            self.qvp=(qv*self.topomask[None,:,:,None]-mqv)
        whereinf=isinf(self.qvp);
        self.qvp[whereinf] = nan;
        nanqvp=self.qvp[:,:,:,:]
        wherefin = isfinite(nanqvp);
        nanqvp[wherefin==False] = 0.0;
        # integrated/min/max quantities
        self.buoymax=nanmax(nanbuoy,axis=1)
        self.buoymin=nanmin(nanbuoy,axis=1)
        self.posbuoypath=sum(nanbuoy*(nanbuoy>0.0)*self.dh,axis=1,dtype=float64)
        self.negbuoypath=sum(nanbuoy*(nanbuoy<0.0)*self.dh,axis=1,dtype=float64)
        self.poswvapanompath=sum(nanqvp*(nanqvp>0.0)*self.dh,axis=1,dtype=float64)
        self.negwvapanompath=sum(nanqvp*(nanqvp<0.0)*self.dh,axis=1,dtype=float64)
        super(ncheighthelper,self).update()
    def make_topo(self):
        # make masks for the topography
        alt=self.gdim('altitude')
        lalt=len(alt)
        tim=self.gdim('time')
        lent=len(tim)
        if(self.geom=='xz'):
            self.topomask=zeros((lalt,shape(hlower)[1]),int)
            self.hhsmall=zeros((lent,lalt,shape(hlower)[1]))
            for k in xrange(lalt):
                for j in xrange(shape(hlower)[1]):
                    self.topomask[k,j]=1.0*(alt[k]>hlower[0,j])              
                    self.hhsmall[:,k,j]=alt[k]*(alt[k]>hlower[0,j])
        elif(self.geom=='yz'):
            self.topomask=zeros((lalt,shape(hlower)[0]),int)
            self.hhsmall=zeros((lent,lalt,shape(hlower)[0]))
            for k in xrange(lalt):
                for j in xrange(shape(hlower)[0]):
                    self.topomask[k,j]=1.0*(alt[k]>hlower[j,0])
                    self.hhsmall[:,k,j]=alt[k]*(alt[k]>hlower[j,0])
        self.topomasknan=zeros(shape(self.topomask))
        self.topomasknan[:]=self.topomask[:]
        self.topomasknan[self.topomask==0]=np.nan
        # calculate height difference between cells, make sure height
        # below the topography is not counted
        if(self.geom=='xz'):
            h1=(self.hhsmall[:,0,:])[:,None,None,:]            
            h2=0.5*(self.hhsmall[:,1:,:]+self.hhsmall[:,:-1,:])[:,:,None,:]
            h3=(self.hhsmall[:,-1,:])[:,None,None,:]
            h=np.concatenate((h1,h2,h3),axis=1)
            replace_h_hlower(h,hlower,shape(h)[0],shape(h)[1],shape(h)[2],shape(h)[3])
        elif(self.geom=='yz'):
            h1=(self.hhsmall[:,0,:])[:,None,None,:]            
            h2=0.5*(self.hhsmall[:,1:,:]+self.hhsmall[:,:-1,:])[:,:,:,None]
            h3=(self.hhsmall[:,-1,:])[:,None,:,None]
            h=np.concatenate((h1,h2,h3),axis=1)
            replace_h_hlower(h,hlower,shape(h)[0],shape(h)[1],shape(h)[2],shape(h)[3])
        self.dh=h[:,1:,:,:]-h[:,:-1,:,:]        
                            
class ncobject(object,get_variable_class):
    # class for writing to netcdf
    def __init__(self,outfile):
        self.data=[] # links to input data
        self.outvars={}
        self.ncoutname=outdir+outfile
        try:
            os.remove(self.ncoutname)
        except:
            pass
        self.outfile=Dataset(self.ncoutname,'w',format='NETCDF4',zlib=lzlib)
        self.outfile.createDimension('time',0)
        self.outfile.createVariable('time', 'f8', ('time',),zlib=lzlib)
        self.dirvars=[]
        self.dervars=[]
        self.dervarsunits=[]
        self.outfile.close()
        self.tstep=0     
    # get a 2d variable from file
    def gv_2d(self,key):
        try:
            return self.data.variables[(key)][:,nboundlines:-nboundlines,nboundlines:-nboundlines]
        except:
            return([nan])
    # get the units
    def gu(self,key):
        try:
            return self.data.variables[(key)].units
        except:
            return([''])
    # functions below are defined in derived classes when needed
    def makedims(self):
        pass
    def app_dirvars(self):
        pass
    def app_dervars(self):
        pass
    def calc_masks(self):
        pass
    # LAYOUT OF POSTPROCESSING TIME STEP
    def appvars(self):
        self.calc_masks()
        self.app_dirvars()
        self.app_dervars()
    def opener(self,data):       
        self.data=data
        self.varkeys=self.data.variables.keys()  
        self.outfile=Dataset(self.ncoutname,'a',format='NETCDF4',zlib=lzlib)
        if(self.tstep==0):
            self.set_dims()
        self.bt=len(self.outfile.variables['time'])
        self.et=len(self.outfile.variables['time'])+len(self.data.variables['time'][:])
        self.outfile.variables['time'][self.bt:self.et]=self.data.variables['time'][:]
        print(self.outfile.variables['time'][self.bt:self.et])
    def closer(self):
        self.tstep=self.tstep+1
        self.outfile.close()
    def app_tstep(self,data):       
        self.opener(data)
        self.appvars()
        self.closer()
    # functions to initialize dimensions in output
    def set_dims(self):
        pass
    def init_dim(self,dimname,dimvalues):
        self.outfile.createDimension(dimname,len(dimvalues))
        var=self.outfile.createVariable(dimname, 'f8', (dimname,),zlib=lzlib)
        var[:]=dimvalues
    def init_dimz(self):
        self.init_dim('z',self.gdim('altitude'))
    def init_dimx(self):
        self.xs=[round(x) for x in 1000*(self.gdim('rlon')[nboundlines:-nboundlines]-nanmin(self.gdim('rlon')[nboundlines:-nboundlines]))/0.00899289]
        self.init_dim('x',self.xs)
    def init_dimy(self):
        self.ys=[round(x) for x in 1000*(self.gdim('rlat')[nboundlines:-nboundlines]-nanmin(self.gdim('rlat')[nboundlines:-nboundlines]))/0.00899289]
        self.init_dim('y',self.ys)
    def init_dimchannels(self):
        self.init_dim('nsynmsg',self.gdim('nsynmsg'))
    # functions to initialize variables in output (dir=DIRECT, der=DERIVED)
    def try_init_dirvar(self,var,dims,mask=''):
        if(self.tstep==0):
            so=self.outfile.createVariable(var+mask, 'f4', dims,zlib=lzlib)
            so.missing_value = nan
            so.long_name=self.data.variables[var].long_name
            so.units=self.gu(var)
    def try_init_der_var(self,var,dims,mask='',vtype=''):
        if(self.tstep==0):
            so=self.outfile.createVariable(var+mask, 'f4', dims,zlib=lzlib)
            so.missing_value = nan
            so.long_name=var
            if(vtype=='maskfrac'):
                so.units='-'
            elif(vtype=='mf'):
                so.units='kg m s-1'
            elif(vtype=='spec'):
                try:
                    so.units=self.dervarsunits[var]+' m-1'
                except:
                    so.units=self.gu(var)+' m-1'
            else:                
                so.units=self.dervarsunits[var]
    # make all sampled variables
    def make_sampvars(self,insampvars):
        self.sampvars={}
        for i in self.dirvars+self.dervars:
            if i in insampvars:             
                self.sampvars[i]=True
            else:
                self.sampvars[i]=False
    def make_specvars(self,inspecvars):
        self.specvars={}
        for i in self.dirvars+self.dervars:
            if i in inspecvars:             
                self.specvars[i]=True
            else:
                self.specvars[i]=False
    # most low level way to write a field
    def put_var(self,var,field):
        self.outfile.variables[var][self.bt:self.et]=field

# general class xy cross section variables
class statgroup_spectra(ncobject):    
    def __init__(self,geom,outfile):
        super(statgroup_spectra,self).__init__(outfile)
        self.geom=geom
        self.initiated=False
                                        
# general class xy cross section variables
class statgroupint(ncobject):    
    def __init__(self,outfile):
        super(statgroupint,self).__init__(outfile)
    def app_dirvars(self):
        for var in self.dirvars:
            if var in self.varkeys:
               self.try_init_dirvar(var,('time','y','x',))
               self.put_var(var,self.gv_2d(var)) 
    def set_dims(self):
        self.init_dimx()
        self.init_dimy()

class statgroupintlevel(statgroupint):
    def robust_minimum_finder(self,field):
        isminimum=nan*zeros(shape(field))
        self.robust_minimum_weaver(isminimum,field,shape(field)[0],shape(field)[1],shape(field)[2],shape(field)[3])
        return isminimum
    def robust_minimum_weaver(self,isminimum,f,tmax,kmax,jmax,imax):
        code = """
        int i, j, k, t;
        for (t=0; t<tmax; t++) {
            for (k=3; k<kmax-3; k++) {
               for (j=0; j<jmax; j++) {
                  for (i=0; i<imax; i++) {
                     if(f(t,k,j,i)>f(t,k-1,j,i) and f(t,k,j,i)>f(t,k-2,j,i) and f(t,k,j,i)>f(t,k-3,j,i) and f(t,k,j,i)>f(t,k+1,j,i) and f(t,k,j,i)>f(t,k+2,j,i) and f(t,k,j,i)>f(t,k+3,j,i)) { 
                         isminimum(t,k,j,i)=1;
                     }
                  }     
               }
            }
        }
        """
        weave.inline(code,['isminimum','f','tmax','kmax','jmax','imax'],type_converters=converters.blitz,compiler='gcc')  
    def __init__(self):
        super(statgroupintlevel,self).__init__('intlv.'+marker+'.nc')
        self.dirvars=dirintvarslevel
        self.dervars=derintvarslevel
        self.dervarsunits=derintvarsunitslevel
        self.helper=levelhelper       
    def app_dervars(self):
        for var in self.dervars:
            self.try_init_der_var(var,('time','y','x',))
        # calculate and put the actual derived variables
        self.put_var('VWP',self.helper.vwp)
        self.put_var('CWP',self.helper.cwp)
        self.put_var('RWP',self.helper.rwp)
        self.put_var('GWP',self.helper.gwp)
        self.put_var('SWP',self.helper.swp)
        self.put_var('IWP',self.helper.iwp)
        self.put_var('WMAX',self.helper.wmax)
        self.put_var('WMIN',self.helper.wmin)
        self.put_var('SSHF',self.helper.sshf)
        self.put_var('SLHF',self.helper.slhf)
        self.put_var('CLDQCIFRAC',self.helper.cldqcifrac)
        self.put_var('CLDQCIW1FRAC',self.helper.cldqciw1frac)
        self.put_var('CLDTOP',nanmax(self.helper.cld*self.helper.hh,axis=1))
        self.put_var('CLDUPDTOP',nanmax(self.helper.cld*(self.helper.wh>1.0)*self.helper.hh,axis=1))
        # find boundary layer using maximum in second derivative
        # algorithm is relatively robust for daytime?
        thetarhograd=(self.helper.thetarhoh[:,2:,:,:]-self.helper.thetarhoh[:,:-2,:,:])/(self.helper.hh[:,2:,:,:]-self.helper.hh[:,:-2,:,:])
        isminimum=self.robust_minimum_finder(thetarhograd)
        hpbl=nanmin(isminimum*self.helper.hh[:,1:-1,:,:],axis=1)
        self.put_var('HPBLTHETA',hpbl)
        self.put_var('DHPBLTHETA',hpbl-hlower)
        thetarhograd2=(thetarhograd[:,1:,:,:]-thetarhograd[:,:-1,:,:])/(self.helper.hh[:,2:-1,:,:]-self.helper.hh[:,1:-2,:,:])
        isminimum=self.robust_minimum_finder(thetarhograd2)
        hpbl=nanmin(isminimum*self.helper.h[:,2:-2,:,:],axis=1)
        self.put_var('HPBLTHETA2',hpbl)
        self.put_var('DHPBLTHETA2',hpbl-hlower)
        # make some space for new variables
        del thetarhograd,thetarhograd2,isminimum
        qtthadv=(self.gv('AQVT_ADV')+self.gv('AQCT_ADV'))
        dh=self.helper.h[:,1:,:,:]-self.helper.h[:,:-1,:,:]
        # qt convergence by advection only
        self.put_var('QTTCONV',sum(qtthadv*(qtthadv>0.)*dh*self.helper.rhoh,axis=1,dtype=float64))
        self.put_var('QTTDIV',sum(qtthadv*(qtthadv<0.)*dh*self.helper.rhoh,axis=1,dtype=float64))
        self.put_var('QTTNET',sum(qtthadv*dh*self.helper.rhoh,axis=1,dtype=float64))
        # d(rho*w)/dz, measure for mass convergence (positive and negative integrals must equal...)
        drhow=-(concatenate(((self.helper.wh*self.helper.rhoh)[:,:,:,:],0.*(self.helper.wh[:,0,:,:]*self.helper.rhoh[:,0,:,:])[:,None,:,:]),axis=1)-concatenate((0.*(self.helper.wh[:,0,:,:]*self.helper.rhoh[:,0,:,:])[:,None,:,:],(self.helper.wh*self.helper.rhoh)[:,:,:,:]),axis=1))
        self.put_var('DRHOWDZPOSINT',sum(drhow*(drhow>0.),axis=1,dtype=float64))
        self.put_var('DRHOWDZNEGINT',sum(drhow*(drhow<0.),axis=1,dtype=float64))        
        self.put_var('LLWIND',self.helper.llwind)        
                
class statgroupintheight(statgroupint):
    def __init__(self):
        super(statgroupintheight,self).__init__('intz.'+marker+'.nc')
        self.dirvars=dirintvarsheight
        self.dervars=derintvarsheight
        self.dervarsunits=derintvarsunitsheight
        self.helper=heighthelper
    def app_dervars(self):
        for var in self.dervars:
            self.try_init_der_var(var,('time','y','x',))
        self.put_var('BUOYMAX',self.helper.buoymax)
        self.put_var('BUOYMIN',self.helper.buoymin)
        self.put_var('POSBUOYPATH',self.helper.posbuoypath)
        self.put_var('NEGBUOYPATH',self.helper.negbuoypath)
        self.put_var('BUOYPATH',self.helper.posbuoypath+self.helper.negbuoypath)
        self.put_var('POSWVAPANOMPATH',self.helper.poswvapanompath)
        self.put_var('NEGWVAPANOMPATH',self.helper.negwvapanompath)
        self.put_var('WVAPANOMPATH',self.helper.poswvapanompath+self.helper.negwvapanompath)

# general class for domain averaged variables
class statgroup_dom(ncobject):    
    def __init__(self,outfile):
        super(statgroup_dom,self).__init__(outfile)
    def app_dirvars(self):
        for var in self.dirvars:
            if var in self.varkeys:
               self.try_init_dirvar(var,('time',))
               self.put_var(var,mean2d(self.gv_2d(var)))

class statgroup_domlevel(statgroup_dom):
    def __init__(self):
        super(statgroup_domlevel,self).__init__('domlv.'+marker+'.nc')
        self.dirvars=dirdomvarslevel
        self.dervars=derdomvarslevel
        self.dervarsunits=derdomvarsunitslevel
        self.helper=levelhelper
    def app_dervars(self):
        for var in self.dervars:
            self.try_init_der_var(var,('time',))
        # calculate and put the actual derived variables   
        self.put_var('VWP',mean2d(self.helper.vwp))
        self.put_var('CWP',mean2d(self.helper.cwp))
        self.put_var('RWP',mean2d(self.helper.rwp))
        self.put_var('GWP',mean2d(self.helper.gwp))
        self.put_var('SWP',mean2d(self.helper.swp))
        self.put_var('IWP',mean2d(self.helper.iwp))
        self.put_var('WMAX',nanmax(nanmax(self.helper.wmax,axis=2),axis=1))
        self.put_var('WMIN',nanmin(nanmin(self.helper.wmin,axis=2),axis=1))
        self.put_var('SSHF',mean2d(self.helper.sshf))
        self.put_var('SLHF',mean2d(self.helper.slhf))
        self.put_var('CLDQCIFRAC',mean2d(self.helper.cldqcifrac))
        self.put_var('CLDQCIW1FRAC',mean2d(self.helper.cldqciw1frac))
          
class statgroup_domheight(statgroup_dom):
    def __init__(self):
        super(statgroup_domheight,self).__init__('domz.'+marker+'.nc')
        self.dirvars=dirdomvarsheight
        self.dervars=derdomvarsheight
        self.dervarsunits=derdomvarsunitsheight
        self.helper=heighthelper
    def app_dervars(self):
        for var in self.dervars:
            self.try_init_der_var(var,('time',))
        self.put_var('BUOYMAX',nanmax(nanmax(self.helper.buoymax,axis=2),axis=1))
        self.put_var('BUOYMIN',nanmin(nanmin(self.helper.buoymin,axis=2),axis=1))

# general class for domain averaged variables along one dimension
# (hoevmoller diagrams)
class statgroup_hov(ncobject):    
    def __init__(self,geom,outfile):
        super(statgroup_hov,self).__init__(outfile)
        self.geom=geom
    def app_dirvars(self):
        for var in self.dirvars:
            if var in self.varkeys:
               if(self.geom=='xz'):
                   self.try_init_dirvar(var,('time','x',))
               elif(self.geom=='yz'):
                   self.try_init_dirvar(var,('time','y',))
               self.put_var(var,mean_1d(self.gv_2d(var),self.geom))
    def set_dims(self):
        if(self.geom=='xz'):
            self.init_dimx()
        elif(self.geom=='yz'):
            self.init_dimy()

class statgroup_hovlevel(statgroup_hov):
    def __init__(self,geom):
        super(statgroup_hovlevel,self).__init__(geom,'hovlv.'+marker+'.nc')
        self.dirvars=dirhovvarslevel
        self.dervars=derhovvarslevel
        self.dervarsunits=derhovvarsunitslevel
        self.helper=levelhelper
    def app_dervars(self):
        # calculate and put the derived variables   
        if(self.geom=='xz'):
            for var in self.dervars:
                self.try_init_der_var(var,('time','x'))
            self.put_var('WMAX',nanmax(self.helper.wmax,axis=1))
            self.put_var('WMIN',nanmin(self.helper.wmin,axis=1))
        elif(self.geom=='yz'):
            for var in self.dervars:
                self.try_init_der_var(var,('time','y'))
            self.put_var('WMAX',nanmax(self.helper.wmax,axis=2))
            self.put_var('WMIN',nanmin(self.helper.wmin,axis=2))
        self.put_var('VWP',mean_1d(self.helper.vwp,self.geom))
        self.put_var('CWP',mean_1d(self.helper.cwp,self.geom))
        self.put_var('RWP',mean_1d(self.helper.rwp,self.geom))
        self.put_var('GWP',mean_1d(self.helper.gwp,self.geom))
        self.put_var('SWP',mean_1d(self.helper.swp,self.geom))
        self.put_var('IWP',mean_1d(self.helper.iwp,self.geom))
        self.put_var('SSHF',mean_1d(self.helper.sshf,self.geom))
        self.put_var('SLHF',mean_1d(self.helper.slhf,self.geom))
        self.put_var('CLDQCIFRAC',mean_1d(self.helper.cldqcifrac,self.geom))
        self.put_var('CLDQCIW1FRAC',mean_1d(self.helper.cldqciw1frac,self.geom))
        self.put_var('LLWIND',mean_1d(self.helper.llwind,self.geom))
                  
class statgroup_hovheight(statgroup_hov):
    def __init__(self,geom):
        super(statgroup_hovheight,self).__init__(geom,'hovz.'+marker+'.nc')
        self.dirvars=dirhovvarsheight
        self.dervars=derhovvarsheight
        self.dervarsunits=derhovvarsunitsheight
        self.helper=heighthelper
    def app_dervars(self):
        if(self.geom=='xz'):
            for var in self.dervars:
                self.try_init_der_var(var,('time','x'))
            self.put_var('BUOYMAX',nanmax(self.helper.buoymax,axis=1))
            self.put_var('BUOYMIN',nanmin(self.helper.buoymin,axis=1))
        elif(self.geom=='yz'):
            for var in self.dervars:
                self.try_init_der_var(var,('time','y'))
            self.put_var('BUOYMAX',nanmax(self.helper.buoymax,axis=2))
            self.put_var('BUOYMIN',nanmin(self.helper.buoymin,axis=2))
        self.put_var('POSBUOYPATH',mean_1d(self.helper.posbuoypath,self.geom))
        self.put_var('NEGBUOYPATH',mean_1d(self.helper.negbuoypath,self.geom))
        self.put_var('BUOYPATH',mean_1d(self.helper.posbuoypath+self.helper.negbuoypath,self.geom))
        self.put_var('POSWVAPANOMPATH',mean_1d(self.helper.poswvapanompath,self.geom))
        self.put_var('NEGWVAPANOMPATH',mean_1d(self.helper.negwvapanompath,self.geom))
        self.put_var('WVAPANOMPATH',mean_1d(self.helper.poswvapanompath+self.helper.negwvapanompath,self.geom))
                                         
# general class for statistics based on 3d output
# this class is somewhat more extensive
# than the previous examples
class statgroup(ncobject):    
    def __init__(self,geom,outfile):
        super(statgroup,self).__init__(outfile)    
        self.geom=geom
        self.init_specmasks(['upstream','top','downstream','upstreambox','topbox','downstreambox'])
    def app_dirvars(self):
        for var in self.dirvars:
            if var in self.varkeys:
                self.make_var(var)
                if(self.sampvars[var]==True):
                    for mask in self.masks.keys():
                        self.make_var(var,mask=mask)
                self.put_mean_int_mask(var,self.gv(var))
    def app_dervars(self): #different for level and height output
        pass
    def set_dims(self):
        self.init_vdims()
        self.init_hdims()
    def init_vdims(self):
        pass
    def init_hdims(self):
        if(self.geom=='xz'):
            self.init_dimx()
        elif(self.geom=='yz'):
            self.init_dimy()
    def mask_fracs(self):
        for mask in self.masks.keys():
            self.maskfrac[mask]=mean_1d(self.masks[mask].field,self.geom) 
            self.make_dervar(str(mask)+'frac',vtype='maskfrac')
            self.put_mean_int(str(mask)+'frac',1.0*self.masks[mask].field)
    def specmask_fracs(self):
        for mask in self.specmasks.keys():
            self.specmaskfrac[mask]=self.specmasks[mask].field
            self.make_dervar(str(mask)+'frac',vtype='maskfrac')
            if(self.geom=='xz'):
                self.put_mean_int(str(mask)+'frac',1.0*self.specmasks[mask].field[:,:,None])
            elif(self.geom=='yz'):
                self.put_mean_int(str(mask)+'frac',1.0*self.specmasks[mask].field[:,None,:])
    def mask_fracs_2d(self):
        for mask in self.masks.keys():
            self.maskfrac2d[mask]=mean2d(self.masks[mask].field) 
    def init_masks(self,masks):
        self.masks={}
        for maskname in masks:
            self.masks[maskname]=mask()
        self.maskfrac={}
        self.maskfrac2d={}
    def init_specmasks(self,masks):
        self.specmasks={}
        for maskname in masks:
            self.specmasks[maskname]=specmask()
        self.specmaskfrac={}
    def calc_specmask(self):
        # define areas at the top, upsteam, and downstream 
        top = nanmax(hlower)
        bot = nanmin(hlower)
        tres1=bot+0.5*(top-bot)
        tres2=bot+0.01*(top-bot)
        tres3=bot+0.0001*(top-bot)
        if(self.geom=='xz'):
            toploc = (hlower[0,:]).argmax()
            comp=hlower[0,:]          
            locs=(xrange(len(hlower[0,:]))<toploc)         
        elif(self.geom=='yz'):
            toploc = (hlower[:,0]).argmax()
            comp=hlower[:,0]
            locs=(xrange(len(hlower[:,0]))<toploc)
        else:
            return
        upstream=(comp>tres3)*(comp<=tres2)*locs            
        downstream=(comp>tres3)*(comp<=tres2)*(locs==False)
        top=(comp>tres1)
        self.specmasks['upstream'].setfield(upstream[None,:]*(self.hhsmall[:,:]>comp[None,:])*(self.hhsmall[:,:]<(comp[None,:]+500.)))
        self.specmasks['top'].setfield(top[None,:]*(self.hhsmall[:,:]>comp[None,:])*(self.hhsmall[:,:]<(comp[None,:]+500.)))
        self.specmasks['downstream'].setfield(downstream[None,:]*(self.hhsmall[:,:]>comp[None,:])*(self.hhsmall[:,:]<(comp[None,:]+500.)))
        self.specmasks['upstreambox'].setfield(upstream[None,:]*(self.hhsmall[:,:]>comp[None,:])*(self.hhsmall[:,:]<(comp[None,:]+2500.)))
        self.specmasks['topbox'].setfield(top[None,:]*(self.hhsmall[:,:]>comp[None,:])*(self.hhsmall[:,:]<(comp[None,:]+2500.)))
        self.specmasks['downstreambox'].setfield(downstream[None,:]*(self.hhsmall[:,:]>comp[None,:])*(self.hhsmall[:,:]<(comp[None,:]+2500.)))
    def make_spec(self,var,field):
        if(shape(hlower)[0]==1):
           return
        for mask in self.specmasks:
            nsegmentsmax=sum(self.specmasks[mask].field[:,:,:]==True)
            if(nsegmentsmax>0):
                dataout=zeros((nsegmentsmax,shape(field)[2]))
                prepare_spectra_xz(dataout,field,self.specmasks[mask].field[:,:,:],shape(field)[0],shape(field)[1],shape(field)[2],shape(field)[3],nsegmentsmax)
                # calculate distance between grid points
                dy=1000*((nanmax(self.gdim('rlat')))-nanmin(self.gdim('rlat')))/(0.00899289*(len(self.gdim('rlat'))-1))
                (p,wavenr)=asp.spectrum_peri(dataout, Fs=1/dy, pad=False, smooth=False,rmzf=True,scale_by_freq=True)
                if(self.spectra.initiated==False):
                    self.spectra.init_dim('wavenr',wavenr)
                    self.spectra.initiated=True
                self.spectra.try_init_der_var(var,('time','wavenr'),mask=mask,vtype='spec')
                self.spectra.put_var(var+mask,p)
                                                                                                                                                  
# class specific to level statistics
# this includes means in 1d and 2d, cross-sections, conditionally sampled variables
# interpolated output to height levels (for convenience/quick data
# exploration)
class statgroup_level(statgroup):
    # initiatize cross sections, interpolation, spectra
    def __init__(self,geom):
        super(statgroup_level,self).__init__(geom,geom+'lv.'+marker+'.nc')
        self.dirvars=dirvarslevel
        self.dervars=dervarslevel
        self.dervarsunits=dervarsunitslevel
        self.helper=levelhelper
        self.spectra=statgroup_spectra(geom,'spectralv.'+marker+'.nc')
        self.spectra.dervarsunits=dervarsunitslevel    
        self.interp1d=statgroup_heightprof('interp1d.'+marker+'.nc')
        self.interp1d.dervarsunits=dervarsunitslevel    
        self.interp=statgroup_interp(self.geom,'interp.'+marker+'.nc')
        self.interp.dervarsunits=dervarsunitslevel
        self.cross=crossgroup_level(self.geom,'crosslv.'+marker+'.nc')
        self.cross.dervarsunits=dervarsunitslevel
        self.cross.interp=crossgroup_interp(self.geom,'crossinterp.'+marker+'.nc')
        self.cross.interp.dervarsunits=dervarsunitslevel
        self.init_masks(['cld','upd','cldupd','cldupdw1'])
        self.make_sampvars(sampvars)
        self.make_specvars(specvars)
    def calc_masks(self):
        # this is where the actual masks are calculated
        self.masks['cld'].setfield(self.helper.cld)
        self.masks['upd'].setfield((self.helper.wh>0.0))
        self.masks['cldupd'].setfield((self.helper.wh>0.0)*self.helper.cld)
        self.masks['cldupdw1'].setfield((self.helper.wh>1.0)*self.helper.cld)
        self.mask_fracs()
        self.specmask_fracs()
    # Adding a time step (overload because of interp,cross and sepctra)
    def app_tstep(self,data):
        self.interp.opener(data)
        self.interp1d.opener(data)
        self.cross.opener(data)
        self.cross.interp.opener(data)
        self.spectra.opener(data)
        if(self.tstep==0):
            self.interp.init_dim('z',interparr)      
            self.interp1d.init_dim('z',interparr)       
            self.cross.interp.init_dim('z',interparr)      
            self.hsmall=mean_1d(self.helper.h,self.geom)       
            self.hhsmall=mean_1d(self.helper.hh,self.geom)   
            self.calc_specmask()
            self.make_topo()
        super(statgroup_level,self).app_tstep(data)
        self.interp.closer()
        self.interp1d.closer()
        self.cross.closer()
        self.cross.interp.closer()
        self.spectra.closer()
    def make_var(self,var,mask=''):
        if(shape(self.gv(var))[1]==shape(self.gdim('level1'))[0] and mask==''):
            zlev='levelf'
        else:
            zlev='levelh'
        if(self.geom=='xz'):
            self.try_init_dirvar(var,('time',zlev,'x'),mask=mask)
            self.interp.try_init_dirvar(var,('time','z','x'),mask=mask)
            if(mask==''):
                self.cross.try_init_dirvar(var,('time',zlev,'x'),mask=mask)                       
                self.cross.interp.try_init_dirvar(var,('time','z','x'),mask=mask)
        elif(self.geom=='yz'):
            self.try_init_dirvar(var,('time',zlev,'y'),mask=mask)
            self.interp.try_init_dirvar(var,('time','z','y'),mask=mask)
            if(mask==''):
               self.cross.try_init_dirvar(var,('time',zlev,'y'),mask=mask)                       
               self.cross.interp.try_init_dirvar(var,('time','z','y'),mask=mask)
        self.interp1d.try_init_dirvar(var,('time','z'),mask=mask)
    def make_dervar(self,var,mask='',zlev='levelh',vtype=''):
        if(self.geom=='xz'):
            self.try_init_der_var(var,('time',zlev,'x'),mask=mask,vtype=vtype)
            self.interp.try_init_der_var(var,('time','z','x'),mask=mask,vtype=vtype)
            if(mask=='' and vtype==''):
                self.cross.try_init_der_var(var,('time',zlev,'x'),mask=mask,vtype=vtype)
                self.cross.interp.try_init_der_var(var,('time','z','x'),mask=mask,vtype=vtype)
        elif(self.geom=='yz'):
            self.try_init_der_var(var,('time',zlev,'y'),mask=mask,vtype=vtype)
            self.interp.try_init_der_var(var,('time','z','y'),mask=mask,vtype=vtype)
            if(mask=='' and vtype==''):
               self.cross.try_init_der_var(var,('time',zlev,'y'),mask=mask,vtype=vtype)
               self.cross.interp.try_init_der_var(var,('time','z','y'),mask=mask,vtype=vtype)
        self.interp1d.try_init_der_var(var,('time','z'),mask=mask,vtype=vtype)
    def app_dervars(self):
        # append the derived variables
        for var in self.dervars:
             self.make_dervar(var)
             if(self.sampvars[var]==True):
                 for mask in self.masks.keys():
                     self.make_dervar(var,mask=mask)
        dw=deviation_1d(self.gv('W'),self.geom)
        dwplus=dw[:,1:,:,:]  
        dwmin=dw[:,:-1,:,:]
        dwh=0.5*dwplus+dwmin
        wvar=0.5*(dwplus*dwplus+dwmin*dwmin)
        del dwplus,dwmin #decrease mem footprint explicitly
        du=deviation_1d(self.gv('U'),self.geom)
        dv=deviation_1d(self.gv('V'),self.geom)
        if 'cosmo_tke' in marker:
            tkeh=0.5*(self.gv('TKE')[:,1:,:,:]+self.gv('TKE')[:,:-1,:,:])
        else:
            tkeh=0.5*((2.0*self.gv('TKE')[:,1:,:,:])**0.5+(2.0*self.gv('TKE')[:,:-1,:,:])**0.5)        
        self.put_mean_int_mask('TKE',tkeh)         
        self.put_mean_int_mask('LVTTKE',tkeh+0.5*(du*du+dv*dv)+0.5*wvar) # Total TKE wrt 1D mean (hence LV) circulation (includes resolved and sgs)
        self.put_mean_int_mask('LVW2RES',wvar) # vertical component of this, both with and without sgs
        self.put_mean_int_mask('LVW2TOTLES',wvar+(2./3.)*tkeh)   
        del tkeh,wvar
        iexnf=(self.gv('P')/pref)**(-rd/cpd) # inverse exner function      
        theta=self.gv('T')*iexnf
        thl=theta-(rlv/cpd)*self.gv('QC')*iexnf #liquid water potential temperature
        self.put_mean_int_mask('THETA',theta)
        self.put_mean_int_mask('THL',thl)
        dtheta=deviation_1d(theta,self.geom)
        del theta,thl
        self.put_mean_int_mask('QCTOT',self.gv('CLW_CON')*self.gv('CLC_CON')+self.gv('QC')) # qc including (convective!) sgs contribution
        self.put_mean_int_mask('QT',self.helper.qt)
        self.put_mean_int_mask('RHO',self.helper.rhoh)
        self.put_mean_int_mask('THETARHO',self.helper.thetarhoh)
        mthetarhoh=mean_1d(self.helper.thetarhoh,self.geom)
        if self.geom=='xz':
            buoy=grav*(self.helper.thetarhoh-mthetarhoh[:,:,None,:])/mthetarhoh[:,:,None,:]
        else:
            buoy=grav*(self.helper.thetarhoh-mthetarhoh[:,:,:,None])/mthetarhoh[:,:,:,None]         
        self.put_mean_int_mask('LVBUOY',buoy) #buoyancy wrt 1d circulation
        self.put_mean_int_mask('LVBUOYWRES',self.helper.rhoh*dwh*buoy) #buoyancy-flux wrt 1d circulation
        del buoy
        dqt=deviation_1d(self.helper.qt,self.geom)
        self.put_mean_int_mask('RHOW',self.helper.rhoh*self.helper.wh)  # mass flux
        self.put_mean_int_mask('HH',self.helper.hh)
        self.put_mean_int_mask('LVRHOUWRES',self.helper.rhoh*dwh*du) # momentum fluxes
        self.put_mean_int_mask('LVRHOVWRES',self.helper.rhoh*dwh*dv)
        self.put_mean_int_mask('LVRHOTHETAWRES',self.helper.rhoh*dwh*dtheta)
        self.put_mean_int_mask('LVRHOQTWRES',self.helper.rhoh*dwh*dqt)
        self.put_mean_int_mask('CLC_RES',self.helper.cld) # resolved cloud cover
        self.put_mean_int_mask('CLC_TOT',self.helper.cld+self.gv('CLC_CON')*(1.0-self.helper.cld)) # only counts sgs clouds when no grid scale clouds
        # total advective tendencies
        if lbud:
            self.put_mean_int_mask('AQVT_ADVTOT',self.gv('AQVT_ZADV')+self.gv('AQVT_ADV'))
            self.put_mean_int_mask('AQTT_ADVTOT',self.gv('AQVT_ZADV')+self.gv('AQVT_ADV')+self.gv('AQCIT_ZADV')+self.gv('AQCIT_ADV'))
        del dtheta,dqt,du,dv,dwh
        # calculate q_sat and RH
        psat = b1 * exp( b2w*(self.gv('T')-b3)/(self.gv('T')-b4w) )
        qsat = (rd/rv)*psat/(self.gv('P')-(1.0-rd/rv)*psat)
        self.put_mean_int_mask('QSAT',qsat)
        self.put_mean_int_mask('RH',self.helper.qt/qsat)
        del psat,qsat
        psati = b1 * exp( b2i*(self.gv('T')-b3)/(self.gv('T')-b4i) )
        qsati = (rd/rv)*psati/(self.gv('P')-(1.0-rd/rv)*psati)
        self.put_mean_int_mask('QSATI',qsati)
        self.put_mean_int_mask('RHI',self.helper.qt/qsati)
        del psati,qsati
        if lbud:        
            for proc in proclist:
                self.put_mean_int_mask('QTT_'+proc,self.gv('QVT_'+proc)+self.gv('QCIT_'+proc))
                self.put_mean_int_mask('AQTT_'+proc,self.gv('AQVT_'+proc)+self.gv('AQCIT_'+proc))
        # resolved mass-fluxes using mask criteria
        for mask in self.masks.keys():
            self.make_dervar(str(mask)+'mf',vtype='mf')
            self.put_mean_int(str(mask)+'mf',self.masks[mask].field*self.helper.rhoh*self.helper.wh)
    def make_topo(self):
        lenxory=shape(self.hsmall)[2]
        lent=shape(self.hsmall)[0]
        lenout=len(interparr)
        self.interphsmall=zeros((lent,lenout,lenxory))
        self.interphhsmall=zeros((lent,lenout,lenxory))
        for t in xrange(lent):
            for i in xrange(lenxory):
                self.interphhsmall[t,:,i]=1.0*(interparr[:]>self.hhsmall[t,-1,i])
        for t in xrange(lent):
            for i in xrange(lenxory):
                self.interphsmall[t,:,i]=1.0*(interparr[:]>self.hsmall[t,-1,i])
        self.interphsmall[self.interphsmall==0.0]=np.nan
        self.interphhsmall[self.interphhsmall==0.0]=np.nan
    # higher level put_var routines that do interpolation as well
    def put_mean_int(self,var,field):
        self.put_var_int(var,mean_1d(field,self.geom))
    def put_var_int(self,var,field):
        self.put_var(var,field)
        lenxory=shape(field)[2]
        lenz=shape(field)[1]
        lent=shape(field)[0]
        lenout=len(interparr)
        varinterp=zeros((lent,lenout,lenxory))
        if(lenz==shape(self.hhsmall)[1]):
            int_to_height(varinterp,interparr,field,self.hhsmall,lent,lenz,lenxory,lenout)
            self.interp.put_var(var,varinterp*self.interphhsmall)
            self.interp1d.put_var(var,nansum(varinterp*self.interphhsmall,axis=2)/nansum(self.interphhsmall,axis=2))
        else:
            int_to_height(varinterp,interparr,field,self.hsmall,lent,lenz,lenxory,lenout)
            self.interp.put_var(var,varinterp*self.interphsmall)
            self.interp1d.put_var(var,nansum(varinterp*self.interphsmall,axis=2)/nansum(self.interphsmall,axis=2))
    # higher level put_var routines that do interpolation and masks
    def put_var_int_mask(self,var,field,mask):
        self.put_var(var,field/self.maskfrac[mask])
        lenxory=shape(field)[2]
        lenz=shape(field)[1]
        lent=shape(field)[0]
        lenout=len(interparr)
        fminterp=zeros((lent,lenout,lenxory))
        mfinterp=zeros((lent,lenout,lenxory))
        varinterp=zeros((lent,lenout,lenxory))
        int_to_height(varinterp,interparr,field/self.maskfrac[mask],self.hhsmall,lent,lenz,lenxory,lenout)
        int_to_height(fminterp,interparr,field,self.hhsmall,lent,lenz,lenxory,lenout)
        int_to_height(mfinterp,interparr,self.maskfrac[mask],self.hhsmall,lent,lenz,lenxory,lenout)
        self.interp.put_var(var,varinterp*self.interphhsmall)
        self.interp1d.put_var(var,nansum(fminterp,axis=2)/nansum(mfinterp,axis=2))
    def put_mean_int_mask(self,var,field):
        self.put_var_int(var,mean_1d(field,self.geom))
        self.crossput_mean_int(var,extract_1d(field,self.geom))
        if(self.sampvars[var]==True):
            for mask in self.masks.keys():
                if(shape(field)[1]==shape(self.hhsmall)[1]):
                    self.put_var_int_mask(var+mask,mean_1d(field*self.masks[mask].field,self.geom),mask)
                else:
                    f=0.5*(field[:,1:,:,:]+field[:,:-1,:,:])
                    self.put_var_int_mask(var+mask,mean_1d(f*self.masks[mask].field,self.geom),mask)
        if(self.specvars[var]==True):
            if(shape(field)[1]==shape(self.hhsmall)[1]):
                self.make_spec(var,field)
            else:
                self.make_spec(var,field[:,1:,:,:])
    # higher level put_var routines that do interpolation as well, for cross sections
    # could probably be merged with put_mean_int routine
    def crossput_mean_int(self,var,field):
        pass
        self.cross.put_var(var,field)
        lenxory=shape(field)[2]
        lenz=shape(field)[1]
        lent=shape(field)[0]
        lenout=len(interparr)
        varinterp=zeros((lent,lenout,lenxory))
        if(lenz==shape(self.hhsmall)[1]):
            int_to_height(varinterp,interparr,field,self.hhsmall,lent,lenz,lenxory,lenout)
            self.cross.interp.put_var(var,varinterp*self.interphhsmall)
        else:
            int_to_height(varinterp,interparr,field,self.hsmall,lent,lenz,lenxory,lenout)
            self.cross.interp.put_var(var,varinterp*self.interphsmall)
    def init_vdims(self):  
        self.init_dim('levelf',self.gdim('level1'))
        self.init_dim('levelh',self.gdim('level'))
            
# class specific to height statistics    
class statgroup_height(statgroup):
    def __init__(self,geom):
        super(statgroup_height,self).__init__(geom,geom+'z.'+marker+'.nc')
        self.dirvars=dirvarsheight
        self.dervars=dervarsheight
        self.dervarsunits=dervarsunitsheight
        self.helper=heighthelper
        self.spectra=statgroup_spectra(geom,'spectraz.'+marker+'.nc')
        self.spectra.dervarsunits=dervarsunitsheight    
        self.v1d=statgroup_heightprof('prof1d.'+marker+'.nc')
        self.v1d.dervarsunits=dervarsunitsheight
        self.cross=crossgroup_height(self.geom,'crossz.'+marker+'.nc')
        self.cross.dervarsunits=dervarsunitsheight
        self.init_masks(['cld','cldcr','upd','cldupd','cldupdw1'])
        self.make_sampvars(sampvars)
        self.make_specvars(specvars)
    def calc_masks(self):
        w=self.gv('W')
        # this is where the actual masks are calculated
        self.masks['cld'].setfield(self.helper.cld)
        self.masks['cldcr'].setfield(self.helper.cld*(self.deviation_2d(self.helper.rhoh)<0.0))
        self.masks['upd'].setfield((w>0.0))
        self.masks['cldupd'].setfield(self.helper.cld*(w>0.0))
        self.masks['cldupdw1'].setfield(self.helper.cld*(w>1.0))
        self.mask_fracs()
        self.specmask_fracs()
        self.mask_fracs_2d()
    def app_tstep(self,data):
        self.v1d.opener(data)
        self.cross.opener(data)
        self.spectra.opener(data)
        if(self.tstep==0):
            self.v1d.init_dimz()
            self.hhsmall=self.helper.hhsmall
            self.calc_specmask()        
        super(statgroup_height,self).app_tstep(data)    
        self.v1d.closer()
        self.cross.closer()
        self.spectra.closer()
    def make_var(self,var,mask=''):
        if(self.geom=='xz'):
            self.try_init_dirvar(var,('time','z','x'),mask=mask)
            if(mask==''):
                self.cross.try_init_dirvar(var,('time','z','x'),mask=mask)
        elif(self.geom=='yz'):
            self.try_init_dirvar(var,('time','z','y'),mask=mask)
            if(mask==''):
                self.cross.try_init_dirvar(var,('time','z','y'),mask=mask)
        self.v1d.try_init_dirvar(var,('time','z'),mask=mask)
    def make_dervar(self,var,mask='',vtype=''):
        if(self.geom=='xz'):
            self.try_init_der_var(var,('time','z','x'),mask=mask,vtype=vtype)
            if(mask=='' and vtype==''):
                self.cross.try_init_der_var(var,('time','z','x'),mask=mask,vtype=vtype)
        elif(self.geom=='yz'):
            self.try_init_der_var(var,('time','z','y'),mask=mask,vtype=vtype)
            if(mask=='' and vtype==''):
                self.cross.try_init_der_var(var,('time','z','y'),mask=mask,vtype=vtype)
        self.v1d.try_init_der_var(var,('time','z'),mask=mask,vtype=vtype)
    def app_dervars(self):
        for var in self.dervars:
             self.make_dervar(var)
             if(self.sampvars[var]==True):
                 for mask in self.masks.keys():
                     self.make_dervar(var,mask=mask)
        w=self.gv('W')
        du=self.deviation_2d(self.gv('U'))
        dv=self.deviation_2d(self.gv('V'))
        dw=self.deviation_2d(w)
        dqt=self.deviation_2d(self.helper.qt)
        iexnf=(self.gv('P')/pref)**(-rd/cpd) #inverse exner function       
        theta=self.gv('T')*iexnf
        thl=theta-(rlv/cpd)*self.gv('QC')*iexnf
        dthl=self.deviation_2d(thl)
        tke=self.gv('TKE')
        ttke=tke+0.5*(du*du+dv*dv+dw*dw) # includes mean circulation
        self.put_mean_int_mask('BUOY',self.helper.buoy)
        self.put_mean_int_mask('TKE',tke) 
        self.put_mean_int_mask('TTKE',ttke) 
        self.put_mean_int_mask('RHO',self.helper.rhoh)
        # variances, * preferred over **2 for computational reasons
        self.put_mean_int_mask('QTP2',dqt*dqt)
        self.put_mean_int_mask('THLP2',dthl*dthl)
        self.put_mean_int_mask('BUOYP2',self.helper.buoy*self.helper.buoy)
        self.put_mean_int_mask('QTP',dqt)
        self.put_mean_int_mask('RHOWBUOY',self.helper.rhoh*w*self.helper.buoy)
        self.put_mean_int_mask('W2',w*w)
        for mask in self.masks.keys():
            self.make_dervar(str(mask)+'mf',vtype='mf')
            self.put_mean_int(str(mask)+'mf',self.masks[mask].field*self.helper.rhoh*self.gv('W'))
    def deviation_2d(self,field):
         if(self.geom=='xz'):
             mfield=mean2d(field*self.helper.topomask[None,:,None,:])[:,:,None,None]/mean2d(self.helper.topomask[None,:,None,:])[:,:,None,None]
             pfield=(field*self.helper.topomask[None,:,None,:]-mfield)        
         elif(self.geom=='yz'):
             mqv=mean2d(field*self.helper.topomask[None,:,:,None])[:,:,None,None]/mean2d(self.helper.topomask[None,:,:,None])[:,:,None,None]
             pfield=(field*self.helper.topomask[None,:,:,None]-mfield)
         return pfield
    def put_mean_int(self,var,field):
        if(self.geom=='xz'):      
            self.put_var(var,mean_1d(field*self.helper.topomask[None,:,None,:],self.geom)*self.helper.topomasknan[None,:,:])
            self.v1d.put_var(var,mean2d(field*self.helper.topomask[None,:,None,:]))
        elif(self.geom=='yz'):
            self.put_var(var,mean_1d(field*self.helper.topomask[None,:,:,None],self.geom)*self.helper.topomasknan[None,:,:])
            self.v1d.put_var(var,mean2d(field*self.helper.topomask[None,:,:,None]))        
    def put_mean_int_mask(self,var,field):
        self.put_var(var,mean_1d(field,self.geom)*self.helper.topomasknan)
        self.cross.put_var(var,extract_1d(field,self.geom)*self.helper.topomasknan)
        self.v1d.put_var(var,mean2d(field))
        if(self.sampvars[var]==True):
            for mask in self.masks.keys():
                self.put_var(var+mask,mean_1d(field*self.masks[mask].field,self.geom)/self.maskfrac[mask])
                self.v1d.put_var(var+mask,mean2d(field*self.masks[mask].field)/self.maskfrac2d[mask])
        if(self.specvars[var]==True):
            if(shape(field)[1]==shape(self.hhsmall)[1]):
                self.make_spec(var,field)
            else:
                self.make_spec(var,field[:,1:,:,:])
    def init_vdims(self):
        self.init_dimz()

# general class for 2d profiles,cross sections, interpolated profiles
# actual processing is part of the statgroup routines
# (for computational reasons)
class statgroup_heightprof(ncobject):    
    pass

class statgroup_interp(statgroup):    
    pass

class crossgroup_level(statgroup):    
    def init_vdims(self):  
        self.init_dim('levelf',self.gdim('level1'))
        self.init_dim('levelh',self.gdim('level'))

class crossgroup_height(statgroup):    
    def init_vdims(self):
        self.init_dimz()
        
class crossgroup_interp(statgroup):    
    pass
    
def detect_geometry():               
    # detect if the geometry of the case is 2-dimensional in any of the directions
    # using HHL[0] of the first .nc file of the levels type
    # hlower will also be used to do the topographic filtering later on
    # also determine how interpolation is done
    global geom,hlower,interparr
    thresh=0.01
    print filelist['level'][0]
    try:
        hlower=var_from_file(Dataset(filelist['level'][0],'r'),'HHL')[0,-1,:,:]
    except:
        print 'could not detect geometry'
        raise
    if(var(hlower[:,0])<thresh):
        geom='xz'
    elif(var(hlower[0,:])<thresh):
        geom='yz'
    else:
        geom='2d'      
    print('the geometry of this experiment is detected as '+geom)   
    # special case for 2d experiments
    if(shape(hlower)[0]==1):
        geom='xz'
    # derive heights where interpolation needed
    hmax=nanmax(var_from_file(Dataset(filelist['level'][0],'r'),'HHL'))
    # directly take input heights if no interpolation needed
    if(var(hlower[:,:])<thresh):
        interparr=var_from_file(Dataset(filelist['level'][0],'r'),'HHL')[0,::-1,0,0]
    elif(hmax>12000.):
        interparr=arange(0.,hmax,100.)
    elif(hmax>6000.):
        interparr=arange(0.,hmax,50.)
    elif(hmax>2000.):
        interparr=arange(0.,hmax,20.)
    else:
        interparr=arange(0.,hmax,10.)
                
# process level based data
def process_levelbased():
    global levelhelper
    levelhelper=nclevelhelper(geom)
    levelout=statgroup_level(geom)
    intlevelout=statgroupintlevel()
    domlevelout=statgroup_domlevel()
    hovlevelout=statgroup_hovlevel(geom)
    if(geom == '2d'):
        print('no level-statistics are created')
    for file_to_process in filelist['level']:
        cosmodata=Dataset(file_to_process,'r',format='NETCDF4')
        levelhelper.update(cosmodata)
        if(geom != '2d'):
            levelout.app_tstep(cosmodata)
        intlevelout.app_tstep(cosmodata)
        domlevelout.app_tstep(cosmodata)
        hovlevelout.app_tstep(cosmodata)
        print clock()-start

# process height based data
def process_heightbased():
    global heighthelper
    heighthelper=ncheighthelper(geom)
    heightout=statgroup_height(geom)
    intheightout=statgroupintheight()
    domheightout=statgroup_domheight()
    hovheightout=statgroup_hovheight(geom)
    for file_to_process in filelist['height']:
        cosmodata=Dataset(file_to_process,'r',format='NETCDF4')
        heighthelper.update(cosmodata)
        heightout.app_tstep(cosmodata)
        intheightout.app_tstep(cosmodata)
        domheightout.app_tstep(cosmodata)
        hovheightout.app_tstep(cosmodata)
        print clock()-start
     
# replace missing values for reading in ncview
# and copy to project (storage) directory
def copy_files_to_project():
    outfiles=glob.glob(outdir+'*'+marker+'*.nc')
    outfilescloud=glob.glob(outdir+'clouds/*.'+marker+'.nc')
    for outfile in outfiles+outfilescloud:
        print outfile
        outfiledata=Dataset(outfile,'r+',format='NETCDF4')
        for var in outfiledata.variables:
            if(type(outfiledata.variables[(var)])==np.ma.masked_array):
                vardata=outfiledata.variables[(var)][:]
                vardata=vardata.filled(nan)
                whereinf=isinf(vardata);
                vardata[whereinf]=nan
        outfiledata.close()
        print clock()-start
    for outfile in outfiles:        
        shutil.copy(outfile,projectdir)
    make_tarfile(outdir+'clouds.'+marker+'.tar',glob.glob(outdir+'clouds/*.'+marker+'.nc'))
    shutil.copy(outdir+'clouds.'+marker+'.tar',projectdir)

# update the variables to post-process by level type
def update_variables():
    global dirvarslevel,dirvarsheight
    global dervarslevel,dervarsheight
    global dirintvarslevel,dirintvarsheight
    global derintvarslevel,derintvarsheight
    global dirhovvarslevel,dirhovvarsheight
    global derhovvarslevel,derhovvarsheight
    global derdomvarslevel,derdomvarsheight
    global derintvarsunitslevel,derintvarsunitsheight
    global dervarsunitslevel
    global sampvars
    
    #add budget variables
    if lbud:
        for proc in proclist:
            dirvarslevel+=['TT_'+proc,'ATT_'+proc,'QVT_'+proc,'AQVT_'+proc,'QCIT_'+proc,'AQCIT_'+proc,]
            dervarsunitslevel.update({'AQTT_'+proc:u'kg kg-1 s-1','QTT_'+proc:u'kg kg-1 s-1',})
        sampvars+=[
        'AQTT_MIC',
        'AQTT_TURB',
        'AQTT_ADVTOT',
        'AQTT_HD'
        ]
        dervarsunitslevel.update({'AQVT_ADVTOT':u'kg kg-1 s-1','AQTT_ADVTOT':u'kg kg-1 s-1'})
    
    #direct variables
    dirvarsheight=dirvarslevel

    # add domain mean variables to hovmoeller variables
    dirhovvarslevel+=dirdomvarslevel
    dirhovvarsheight+=dirdomvarsheight
    derhovvarsunitslevel.update(derdomvarsunitslevel)
    derhovvarsunitsheight.update(derdomvarsunitsheight)
             
    # add hovmoeller variables to xy (int) variables
    dirintvarslevel+=dirhovvarslevel
    dirintvarsheight+=dirhovvarsheight
    derintvarsunitslevel.update(derhovvarsunitslevel)
    derintvarsunitsheight.update(derhovvarsunitsheight)
       
    # make list of derived variables separted from keys
    dervarsheight=dervarsunitsheight.keys()
    dervarslevel=dervarsunitslevel.keys()
    derintvarslevel=derintvarsunitslevel.keys()
    derintvarsheight=derintvarsunitsheight.keys()
    derdomvarslevel=derdomvarsunitslevel.keys()
    derdomvarsheight=derdomvarsunitsheight.keys()
    derhovvarslevel=derhovvarsunitslevel.keys()
    derhovvarsheight=derhovvarsunitsheight.keys()

# general class for the crosssections
# take into account areas below topography as nans
# use just one variable per file!
class xycross_onevar(ncobject):    
    def __init__(self,var,outfile):
        super(xycross_onevar,self).__init__(outfile)
        self.dirvars=[var]
    def make_topo(self):
        alt=self.gdim('altitude')
        lalt=len(alt)
        self.topomask=zeros((lalt,shape(hlower)[0],shape(hlower)[1]),int)
        self.topomasknan=zeros((lalt,shape(hlower)[0],shape(hlower)[1]),float)
        for k in xrange(lalt):
            for j in xrange(shape(hlower)[0]):
                for i in xrange(shape(hlower)[1]):
                    self.topomask[k,j,i]=1*(alt[k]>hlower[j,i])
        self.topomasknan=zeros(shape(self.topomask))
        self.topomasknan[:]=self.topomask[:]
        self.topomasknan[self.topomask==0]=np.nan
    def app_dirvars(self):
        if(self.tstep==0):
            self.make_topo()
        for var in self.dirvars:
            if var in self.varkeys:
                self.make_var(var)
                self.put_var(var,self.gv(var)*self.topomasknan)
    def set_dims(self):
        self.init_vdims()
        self.init_hdims()
    def init_vdims(self):
        self.init_dimz()
    def init_hdims(self):
        self.init_dimx()
        self.init_dimy()
    def make_var(self,var):
        self.try_init_dirvar(var,('time','z','y','x'))
                        
def compress_crossxy():
    vardata=Dataset(filelist['crossxy'][0])
    varkeys=vardata.variables.keys()
    compressed_files={}
    for var in varkeys:
        if len(shape(vardata.variables[var]))==4:
           compressed_files[var]=xycross_onevar(var,'crossxy.'+var+'.'+marker+'.nc')
    for file_to_process in filelist['crossxy']:
        cosmodata=Dataset(file_to_process,'r',format='NETCDF4')
        for var in varkeys:
            if len(shape(vardata.variables[var]))==4:
                compressed_files[var].app_tstep(cosmodata)

# general class for compressed files containing satellite output
# here we save with reduced precision and compress a lot
class satvars_compressed(ncobject):    
    def __init__(self,varlist,outfile):
        super(satvars_compressed,self).__init__(outfile)
        self.dirvars=varlist    
    def app_dirvars(self):
        for var in self.dirvars:
            if var in self.varkeys:
                self.make_var(var)
                self.put_var(var,self.gv(var))
    def try_init_dirvar(self,var,dims):
        if(self.tstep==0):
            so=self.outfile.createVariable(var,'f4', dims,zlib=lzlib,least_significant_digit=6)
            so.missing_value = nan
            so.long_name=self.data.variables[var].long_name
            so.units=self.gu(var)
    def set_dims(self):
        self.init_vdims()
        self.init_hdims()
    def init_vdims(self):
        self.init_dimchannels()
    def init_hdims(self):
        self.init_dimx()
        self.init_dimy()
    def make_var(self,var):
        self.try_init_dirvar(var,('time','nsynmsg','y','x'))
                               
def compress_sats():
    if len(filelist['sats'])==0:
        return
    vardata=Dataset(filelist['sats'][0])
    varkeys=vardata.variables.keys()
    vars_to_process=[]
    for var in varkeys:
        if len(shape(vardata.variables[var]))==4:
           vars_to_process.append(var)
    satout=satvars_compressed(vars_to_process,'sats.'+marker+'.nc')
    for file_to_process in filelist['sats']:
        cosmodata=Dataset(file_to_process,'r',format='NETCDF4')
        satout.app_tstep(cosmodata)
                
# general class for compressed files containing cloud stats
# here we save with reduced precision and compress a lot
class cloudvars_compressed(ncobject):    
    def __init__(self,varlist,outfile):
        super(cloudvars_compressed,self).__init__(outfile)
        self.dirvars=varlist    
    def app_dirvars(self):
        for var in self.dirvars:
            if var in self.varkeys:
                self.make_var(var)
                self.put_var(var,self.gv(var))
        self.add_hsurf()
    def try_init_dirvar(self,var,dims):
        if(self.tstep==0):
            so=self.outfile.createVariable(var,'f4', dims,zlib=lzlib,least_significant_digit=6)
            so.missing_value = nan
            so.long_name=self.data.variables[var].long_name
            so.units=self.gu(var)
    def set_dims(self):
        self.init_vdims()
        self.init_hdims()
    def init_vdims(self):
        self.init_dimz()
    def init_hdims(self):
        self.init_dimx()
        self.init_dimy()
    def make_var(self,var):
        self.try_init_dirvar(var,('time','z','y','x'))
    def add_hsurf(self):
        self.try_init_dirvar('HSURF',('time','y','x'))
        self.put_var('HSURF',self.gv_2d('HSURF'))
        
def compress_clouds():
    vardata=Dataset(filelist['crossxy'][0])
    varkeys=vardata.variables.keys()
    vars_to_process=[]
    for var in varkeys:
        if len(shape(vardata.variables[var]))==4 and (var!='QV'): #do not include QV, HSURF and dimensions
           vars_to_process.append(var)
    for file_to_process in filelist['cloud']:
        cloudout=cloudvars_compressed(vars_to_process,'clouds/cloud.'+file_to_process[-13:-5]+'.'+marker+'.nc')
        cosmodata=Dataset(file_to_process,'r',format='NETCDF4')
        cloudout.app_tstep(cosmodata)
                            
########### MAIN PROGRAM ########### 

def runme():
    update_variables()
    mkdir_p(outdir)
    mkdir_p(outdir+'/clouds')
    mkdir_p(projectdir)
    make_filelist()
    detect_geometry()
    process_levelbased()
    process_heightbased()
    if lcross:
        compress_crossxy()
    if lclouds:
        compress_clouds()
    if lsats:
        compress_sats()
    copy_files_to_project()

# ACTUALLY CALLS THE SCRIPT FROM THE COMMAND LINE
# Using: if __name__ == "__main__" 
# makes sure we can import the separate routines

if __name__ == "__main__":
    case=sys.argv[1]
    conf=sys.argv[2]
    exper=sys.argv[3]
    fulldir='/scratch/daint/'+myusername+'/exe/csim/fromtempl/'+case+'/'+conf+'/'+exper+'/output/'
    outdir='/scratch/daint/'+myusername+'/statdump/'+case+'/'+exper+'/'
    projectdir='/project/ch4/'+myusername+'/statdump/'+case+'/'+exper+'/'   
    marker=conf
    runme()
