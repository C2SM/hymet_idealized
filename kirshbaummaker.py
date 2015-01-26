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

from pylab import *
from cosmo_utils import *
from scipy.interpolate import pchip_interpolate 
import os

psurf=101200
pref=100000
r_d=287.05
r_v=461.5
c_pd = 1004.0
grav=9.81
c_pv = 1410
lh_v = 2.501E6
umax = 0.0

Hsurf=0
Hbltop=1000
HmidRH=6250
Hstrat=11500
Hlowerfree=3000

RHsurf=0.7
RHbltop=0.5
RHmid=0.2
RHstrat=0.5

height=22500
RH=zeros(height)
T=zeros(height)
pp=zeros(height)

Tbl=296
Nlowerfree=0.002
Hlowerfree=3000
Nstrat=0.02
grav=9.81
sfratio=0.3
fluxmax=420

pp[0]=psurf
# Boundary layer
for zz in range(Hsurf,Hbltop):
    RH[zz]=RHsurf+(RHbltop-RHsurf)*(zz-Hsurf)/(Hbltop-Hsurf)
    T[zz]=Tbl
    psat=fspw(T[zz])
    qv=RH[zz]*fsqv(psat, pp[zz])
    Tv=fstv(T[zz],qv)
    pp[zz+1]=pp[zz]*exp(-grav/(r_d*Tv))

# Tropospheric and stratospheric RH
qfitRH = polyfit([Hbltop,HmidRH,Hstrat], [RHbltop,RHmid,RHstrat], 2)

for zz in range(Hbltop,Hstrat):
    RH[zz]=qfitRH[0]*zz**2.0+qfitRH[1]*zz+qfitRH[2]
        
for zz in range(Hstrat,height):
    RH[zz]=RHstrat

# Tropospheric and stratospheric T
# Iterate the hydostatic equation a few times
for zz in range(Hbltop,Hlowerfree):
    exn=(pp[zz]/pref)**(r_d/c_pd)
    exnmin=(pp[zz-1]/pref)**(r_d/c_pd)
    thetavmin=Tv/exnmin
    thetav=thetavmin*exp(Nlowerfree**2/grav)
    Tv=thetav*exn
    qtry=qv
    for i in range(5):
        Ttry=Tv/(1 + 0.61*qtry)
        psat=fspw(Ttry)
        qtry=RH[zz]*fsqv(psat, pp[zz])
    T[zz]=Ttry
    psat=fspw(T[zz])
    qv=RH[zz]*fsqv(psat, pp[zz])
    Tv=fstv(T[zz],qv)
    pp[zz+1]=pp[zz]*exp(-grav/(r_d*Tv))

# Use equation from ams glossary for pseudoadiabat
for zz in range(Hlowerfree,Hstrat):
    mrv=qv/(1-qv)
    Ttry=T[zz-1]
    Ttry=T[zz-1]-grav*(1+mrv)*(1+lh_v*mrv/(r_d*Ttry))/(c_pd+c_pv*mrv+(lh_v*lh_v*mrv*(r_d/r_v+mrv))/(r_d*Ttry*Ttry))
    psat=fspw(Ttry)
    qtry=RH[zz]*fsqv(psat, pp[zz])
    Tmean=0.5*(T[zz]+Ttry)
    qmean=0.5*(qv+qtry)
    mrvmean=qmean/(1-qmean)
    for i in range(5):
        Ttry=T[zz-1]-grav*(1+mrvmean)*(1+lh_v*mrvmean/(r_d*Tmean))/(c_pd+c_pv*mrvmean+(lh_v*lh_v*mrvmean*(r_d/r_v+mrvmean))/(r_d*Tmean*Tmean))
        psat=fspw(Ttry)
        qtry=RH[zz]*fsqv(psat, pp[zz])
        Tmean=0.5*(T[zz]+Ttry)
        qmean=0.5*(qv+qtry)
        mrvmean=qmean/(1-qmean)
    T[zz]=Ttry    
    psat=fspw(T[zz])
    qv=RH[zz]*fsqv(psat, pp[zz])
    Tv=fstv(T[zz],qv)
    pp[zz+1]=pp[zz]*exp(-grav/(r_d*Tv))    

for zz in range(Hstrat,height):
    exn=(pp[zz]/pref)**(r_d/c_pd)
    exnmin=(pp[zz-1]/pref)**(r_d/c_pd)
    thetavmin=Tv/exnmin
    thetav=thetavmin*exp(Nstrat**2/grav)
    Tv=thetav*exn
    qtry=qv
    for i in range(5):
        Ttry=Tv/(1 + 0.61*qtry)
        psat=fspw(Ttry)
        qtry=RH[zz]*fsqv(psat, pp[zz])
    T[zz]=Ttry
    psat=fspw(T[zz])
    qv=RH[zz]*fsqv(psat, pp[zz])
    Tv=fstv(T[zz],qv)
    if(zz<height-1):
        pp[zz+1]=pp[zz]*exp(-grav/(r_d*Tv))

zztarg=array(range(0,22500,50))
pt=pchip_interpolate(range(height),T,zztarg)
rh=pchip_interpolate(range(height),RH,zztarg)
uu=umax*tanh(zztarg/2000.)
        
dcosmo = np.zeros((zztarg.size,8))
dcosmo[0,0] = psurf/100.
dcosmo[:,1] = zztarg
dcosmo[:,2] = pt
dcosmo[:,3] = 0.0
dcosmo[:,4] = rh*100.
dcosmo[:,5] = 0.0 
dcosmo[:,6] = uu
dcosmo[:,7] = 270.0 

# create sounding file for COSMO
# ------------------------------
filename = 'output_folder/atmos.input'
print 'writing cosmo format: ' + filename
np.savetxt("snd_cosmo_noheader.input", dcosmo[:,:], fmt='%9.3f')

header_cosmo= "# Input for COSMO\n"
header_cosmo += "# \n"
header_cosmo += "# \n"
header_cosmo += " P [hPa]     Z [m]     T [K]  Dewp [K]    RH [%]  r [g/kg]  WS [m/s]  WD [deg]"

fout = open(filename, 'w')
print >>fout, header_cosmo
fin = open("snd_cosmo_noheader.input")
line = fin.readline()
while line:
    print >>fout, line,
    line = fin.readline()
fin.close()
fout.close()

os.remove('snd_cosmo_noheader.input')

filename = 'output_folder/fluxes.input'
print 'writing cosmo format: ' + filename
header_cosmo= "# Input for COSMO\n"
header_cosmo += "# \n"
header_cosmo += "# \n"
header_cosmo += " time [s]     SHF [W/m**2]     LHF [W/m**2]"
fout = open(filename, 'w')
print >>fout, header_cosmo
for tt in range(0,12*3600+10,10):
    print >>fout, "%5i     %9.4f     %9.4f"%(tt,sfratio*sin(2*pi*tt/(24*3600.))*fluxmax,(1-sfratio)*sin(2*pi*tt/(24*3600.))*fluxmax)
fout.close()
