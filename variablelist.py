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

#processes for budgets
proclist=[
'MIC',
'ADV',
'TURB',
'TOT',
'HD',
'ZADV',
'RLX',
'RELAX',
'LSF',
'CON',
]

dirvarslevel=[
'QV',
'QR',
'QS',
'QI',
'QC',
'QG',
'T',
'U',
'V',
'W',
'TKVH',
'TKVM',
'P',
'EFLUX',
'HFLUX',
'CLC_CON',
'HHL',
'CLW_CON',
'UT_LSF',
'VT_LSF',
'TKESV-QT',
'TKESV-QQ',
'TKESV-TT',
]

dervarsunitslevel={
'THETA':u'K',
'THETARHO':u'K',
'LVTTKE':u'm2 s-2',
'LVW2RES':u'm2 s-2',
'QCTOT':u'kg-1 kg',
'QT':u'kg-1 kg',
'THL':u'K',
'LVW2TOTLES':u'm2 s-2',
'RHO':u'kg',
'RHOW':u'kg m s-1',
'LVBUOY':u'm s-2',
'QSAT':u'kg-1 kg',
'QSATI':u'kg-1 kg',
'RH':u'-',
'RHI':u'-',
'HH':u'm',
'LVRHOUWRES':u'kg m-1 s-2',
'LVRHOVWRES':u'kg m-1 s-2',
'LVRHOTHETAWRES':u'kg K m-2 s-1',
'LVRHOQTWRES':u'kg m-2 s-1',
'LVBUOYWRES':u'm-1 s-1',
'TKE':u'm2 s-2',
'CLC_RES':u'-',
'CLC_TOT':u'-',
}
   
dervarsunitsheight={
'TKE':u'm2 s-2',
'TTKE':u'm2 s-2',
'RHO':u'kg m-3',
'BUOY':u'm s-2',
'QTP2':u'kg-2 kg2',
'THLP2':u'K2',
'BUOYP2':u'K2',
'RHOWBUOY':u'kg m-2 s-1 K',
'W2':u'm2 s-2',
'QTP':u'kg kg-1',
}

### DOM VARS
dirdomvarslevel=[
'RAIN_GSP',
'CAPE_ML',
'CIN_ML',
'CAPE_MU',
'CIN_MU',
'CLCH',
'CLCM',
'CLCL',
'CLCT',
'MFLX_CON',
'HTOP_CON',
'HBAS_CON']

dirdomvarsheight=[]

derdomvarsunitslevel={
'VWP':u'kg m-2',
'CWP':u'kg m-2',
'RWP':u'kg m-2',
'GWP':u'kg m-2',
'SWP':u'kg m-2',
'IWP':u'kg m-2',
'WMAX':u'm s-1',
'WMIN':u'm s-1',
'SSHF':u'W m-2',
'SLHF':u'W m-2',
'CLDQCIFRAC':u'-',
'CLDQCIW1FRAC':u'-',
}

derdomvarsunitsheight={
'BUOYMAX':u'm s-1',
'BUOYMIN':u'm s-1',
}

### ADDITIONAL HOVMOELLER DIAGRAM VARS
dirhovvarslevel=[
'HSURF',
]

dirhovvarsheight=[
]

derhovvarsunitslevel={
'LLWIND':u'm s-1',
}

derhovvarsunitsheight={
'POSBUOYPATH':u'm2 s-2 m-2',
'NEGBUOYPATH':u'm2 s-2 m-2',
'BUOYPATH':u'm2 s-2 m-2',
'POSWVAPANOMPATH':u'kg m-2',
'NEGWVAPANOMPATH':u'kg m-2',
'WVAPANOMPATH':u'kg m-2',
}

### ADDITIONAL INT VARS
dirintvarslevel=[
]

dirintvarsheight=[
]

derintvarsunitslevel={
'CLDTOP':u'm',
'CLDUPDTOP':u'm',
'QTTCONV':u'kg m-2',
'QTTDIV':u'kg m-2',
'QTTNET':u'kg m-2',
'HPBLTHETA':u'm',
'DHPBLTHETA':u'm',
'HPBLTHETA2':u'm',
'DHPBLTHETA2':u'm',
'DRHOWDZPOSINT':u'kg s-1 m-2',
'DRHOWDZNEGINT':u'kg s-1 m-2',
}

derintvarsunitsheight={
}

# list of variables for which to produce sampled statistics
sampvars=[
'QV',
'QR',
'QS',
'QI',
'QC',
'QG',
'T',
'W',
'TKE',
'TKVH',
'TKVM',
'P',
'EFLUX',
'HFLUX',
'THETA',
'LVTTKE',
'LVW2RES',
'QCTOT',
'QT',
'THL',
'LVW2TOTLES',
'RHO',
'RHOW',
'LVBUOY',
'LVRHOUWRES',
'LVRHOVWRES',
'LVRHOTHETAWRES',
'LVRHOQTWRES',
'LVBUOYWRES',
'TTKE',
'RHO',
'BUOY',
'TKESV-QT',
'TKESV-QQ',
'TKESV-TT',
]

# list of variables for which to produce spectra
specvars=[
'U',
'V',
'W',
'TKE',
'QT',
'THL',
'BUOY',
]
