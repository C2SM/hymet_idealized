#!/usr/bin/env python

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
# Copyright (C) 2012 Juerg Schmidli, ETHZ

#
# Library functions for atmos_profile.py
# Juerg Schmidli, IACETH, Dezember 2012
#
# Many COSMO constants can be found in src_setup.f90
#----------------------------------------------------------------

import numpy as np

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.itervalues():
        sp.set_visible(False)


def fspw(T):
    """ Calculate saturation vapor pressure over water in Pa 
    
        Input: temperature in Kelvin
        Output: saturation vapor pressure over water in Pa
    """
    # constants for computing the saturation vapour pressure
    # over water (w) and ice (i) according to Teten's formula
    b1  = 610.78
    b2w = 17.2693882
    b3  = 273.16
    b4w = 35.86

    res = b1*np.exp(b2w*(T-b3)/(T-b4w) )
    return res

def fsqv(zspw, p):
    """ Calculate specific humidity at saturation
    
        zspw:   saturation vapor pressure over water in Pa
        p:      pressure in Pa
        Output: specific humidity (kg/kg)
    """
    # physical constants
    t0      = 273.15
    r_d     = 287.05
    r_v     = 461.51
    rdv     = r_d / r_v
    o_m_rdv = 1.0 - rdv     # one minus rdv 

    res = rdv*zspw/( p - o_m_rdv*zspw )
    return res


def fstv(T, qv):
    """ Calculate virtual temperature in Kelvin
    
        T:      temperature [K]
        qv:     specific humidity (kg/kg)
        Output: virtual temperature [K]
    """
    res = T*(1 + 0.61*qv)
    return res


def p_hydro(psrf, dz, Tv):
    """ Calculate hydrostatic pressure from virtual temperature profile and
        surface pressure
    
        psrf:   surface pressure [Pa]
        dz:     layer thickness [m]
        Tv:     virtual temperature [K]

        Output: hydrostatic pressure [Pa]
    """
    # physical constants
    grav = 9.80665
    r_d  = 287.05

    dz = np.atleast_1d(dz)
    nz = dz.size
    p_hydro = np.zeros(nz+1)
    p_hydro[0] = psrf
    for k in xrange(nz):
        p_hydro[k+1] = p_hydro[k] * np.exp(-grav*dz[k]/(r_d*Tv[k]))
    return p_hydro


debug = False
if debug:
    pref=1000.*100.
    print "T [C]    T[K]   e_sat [hPa]  qv_sat [g/kg]  Tv [K]"
    for T in (273.15, 293.15):
        esat = fspw(T)
        qvs = fsqv(esat, pref)
        Tv = fstv(T, qvs)
        print '%7.2f  %7.2f  %7.2f  %7.2f  %7.2f' % (T-273.15, T, esat/100., qvs*1000., Tv)

    dz = [1000,1000]
    Tv = [290, 270]
    ph = p_hydro(pref, dz, Tv)/100.
    print
    print "hydrostatic pressure"
    print "dz   Tv  p_hydro"
    for (a, b, c) in zip(dz, Tv, ph[1:]):
        print '%7.2f  %7.2f  %7.2f' % (a, b, c)


def wind_prof_linquad(u_low, z_low, u_jet, z_jet, u_top, z_top, zz):
    """ Create a wind profile from given key values. 

Create a wind profile as a combination of a straight line (lower troposphere)
and a parabola (upper troposphere and stratosphere)

functions:
    u1(z) = u_low + B (z-z_low),    B>0
    u2(z) = u_jet - A (z-z_jet)**2, A>0
with    
    u2(z_top) = u_top, thus
    A = (u_jet-u_top)/(z_top-z_jet)**2

conditions:
    (1) u1(z_low) = u_low
    (2) u1(z_m) = u2(z_m) = u_m
    (3) u1'(z_m) = u2'(z_m)

from (1), (2) and (3):
    u_m = u_jet - A(z_m-z_jet)**2
    -2A(u_m-z_jet) = (u_m-u_low)/(z_m-z_low)

results in quadratic equation for x=z_m-z_jet:
   A x^2 + 2A(z_jet-z_low) x + (u_jet-u_low) = 0
with 
    a=A
    b=2A(z_jet-z_low)
    c=u_jet-u_low 
and:
    x1,2 = (-b pm sqrt(b**2-4ac))/(2a)
"""

    # determine A
    A = (u_jet-u_top)/(z_top-z_jet)**2

    # solve quadratic equation
    a = A
    b = 2*A*(z_jet-z_low)
    c = u_jet-u_low
    z_m = z_jet + (-b+np.sqrt(b*b-4*a*c))/(2*a)
    u_m = u_jet - A*(z_m-z_jet)**2
    B = (u_m-u_low)/(z_m-z_low)

    # determine wind profile
    uu = u_jet - A*(zz-z_jet)**2
    idx_lower = np.where(zz < z_m)
    uu[idx_lower] = u_low + B*(zz-z_low)

    return uu


def wind_prof_quad(u_low, z_low, u_jet, z_jet, u_top, z_top, zz):
    """ Create a wind profile from given key values. 

Create a wind profile as a combination of two parabolas 
(after Schlemmer et al, 2011, JAS)

functions:
    u1(z) = u_low + a1 (z-z_low)**2,    a1>0
    u2(z) = u_jet - a2 (z-z_jet)**2,    a2>0
with    
    u1(z_low) = u_low, u2(z_jet) = u_jet, u2(z_top) = u_top, thus
    a2 = (u_jet-u_top)/(z_top-z_jet)**2

matching conditions:
    (1) u1(z_m) = u2(z_m)
    (2) u1'(z_m) = u2'(z_m)
two unknowns: a1, z_m

results in quadratic equation for z_m:
   a z_m**2 + b z_m + c = 0
with 
    a = a2(z_jet - z_low)
    b = u_jet-u_low - a2(z_jet**2 - z_low**2)
    c = z_low*[a2*z_jet*(z_jet-z_low) + u_jet - u_low
and:
    z_m1,2 = (-b pm sqrt(b**2-4ac))/(2a)
and:
    a1 = a2*(z_jet-z_m)/(z_m-z_low)
"""
    # determine a2
    a2 = (u_jet-u_top)/(z_top-z_jet)**2

    # solve quadratic equation
    a = a2*(z_jet - z_low)
    b = u_jet - u_low - a2*(z_jet**2 - z_low**2)
    c = z_low*(a2*z_jet*(z_jet - z_low) + u_jet - u_low)
    z_m = (-b+np.sqrt(b*b-4*a*c))/(2*a)
    a1 = a2*(z_jet-z_m)/(z_m-z_low)
    u_m = u_low + a1*(z_m-z_jet)**2
    #print "a, b, c: ", a, b, c
    #print "z_m, u_m: ", z_m, u_m
    #print "a1: ", a1

    print "matching height and wind speed: %6.2f  %6.2f" % (z_m, u_m)

    # determine wind profile
    uu = u_jet - a2*(zz-z_jet)**2
    idx_lower = np.where(zz < z_m)
    uu[idx_lower] = u_low + a1*(zz-z_low)**2

    idx_lower = np.where(zz < z_low)[0]
    if idx_lower.size > 0:
        uu[idx_lower] = u_low 

    return uu


# zonal wind profile (linear)
#u = u_srf + (u_jet - u_srf)/z_tp*z
#k = (z > z_tp)
#u[k] = u_jet + (u_min - u_jet)/(z_top-z_tp)*z

