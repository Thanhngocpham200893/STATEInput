#!/home/PCoMS/thanhpn/myenv/bin/python3
"""
Read the relaxed cooridinates in GEOMETRY file and prepare the new input file of STATE.
If GEOMETRY is not found (scf calculation), return the old input file. 
Usage: relax_state.py nfinp_old  
"""
#import library here
import numpy as np
from numpy import pi, sin, cos, arccos, sqrt, dot
import sys
import os
##-----end load library


altv = np.zeros((3,3),dtype = float)
b2ang= 0.529177
##read the name of state input from command line

#import library here
import numpy as np
from numpy import pi, sin, cos, arccos, sqrt, dot
import sys
import os
##-----end load library

altv = np.zeros((3,3),dtype = float)
b2ang= 0.529177


##read the name of input files from command line
if len(sys.argv) == 2:
    state_inp     = str(sys.argv[1])
else:
    sys.exit('relax_state.py nfinp_* ')


if len(sys.argv) == 2:
    stateinp = str(sys.argv[1])
else:
    sys.exit('Usage: python3 stateinp2xsf.py nfinp_#.xyz > nfinp_#.xsf')

nfinp = []

with open(stateinp,"r") as inp:
    nfinp_old = inp.readlines()
    for line in nfinp_old:
        if line.startswith('#'):
            continue
        elif line.startswith('\n'):
            continue
        elif len(line.split()) == 0:
            continue
        else:
            nfinp.append(line)

idum1 = int(nfinp[0].split()[0])   
ktyp = int(nfinp[1].split()[2])   
katm = int(nfinp[1].split()[3])

if (idum1 == 1):  
    altv[0,:] =   nfinp[4].split()[0:3] 
    altv[1,:] =   nfinp[5].split()[0:3]  
    altv[2,:] =   nfinp[6].split()[0:3] 
    cps    = np.zeros((katm,3),dtype = float)
    ktyp_l = np.zeros((katm,1),dtype = int)
    start = 9
    for i in range(katm):
        dummy = nfinp[start+i]
        cps[i,0:3] = dummy.split()[0:3]
        ktyp_l[i] = dummy.split()[5]
    z_ia = np.zeros((ktyp,1),dtype = float)
    for i in range(ktyp):
        dummy = nfinp[start+katm+i]
        z_ia[i,0] = dummy.split()[0]
else:
    a     = float(nfinp[3].split()[0])   
    b     = float(nfinp[3].split()[0])   
    c     = float(nfinp[3].split()[0])   
    alpha = float(nfinp[3].split()[0])   
    beta  = float(nfinp[3].split()[0])   
    gamma = float(nfinp[3].split()[0])   
# Handle orthorhombic cells separately to avoid rounding errors, taken from ASE 
    eps = 2 * np.spacing(90.0, dtype=np.float64)  # around 1.4e-14
# alpha
    if abs(abs(alpha) - 90) < eps:
        cos_alpha = 0.0
    else:
        cos_alpha = cos(alpha * pi / 180.0)
# beta
    if abs(abs(beta) - 90) < eps:
        cos_beta = 0.0
    else:
        cos_beta = cos(beta * pi / 180.0)
# gamma
    if abs(gamma - 90) < eps:
        cos_gamma = 0.0
        sin_gamma = 1.0
    elif abs(gamma + 90) < eps:
        cos_gamma = 0.0
        sin_gamma = -1.0
    else:
        cos_gamma = cos(gamma * pi / 180.0)
        sin_gamma = sin(gamma * pi / 180.0)

# Build the cell vectors
    altv[0,:] = a * np.array([1, 0, 0])
    altv[1,:] = b * np.array([cos_gamma, sin_gamma, 0])
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz_sqr = 1. - cx * cx - cy * cy
    assert cz_sqr >= 0
    cz = sqrt(cz_sqr)
    altv[2,:] = c * np.array([cx, cy, cz])

    cps = np.zeros((katm,3),dtype = float)
    ktyp_l = np.zeros((katm,1),dtype = int)
    start = 6
    for i in range(katm):
        dummy = nfinp[start+i]
        cps[i,0:3] = dummy.split()[0:3]
        ktyp_l[i] = dummy.split()[5]
    z_ia = np.zeros((ktyp,1),dtype = float)
    for i in range(ktyp):
        dummy = nfinp[start+katm+i]
        z_ia[i,0] = dummy.split()[0]


###read relaxed coordinate from GEOMETRY
path = "./"
dir_list = os.listdir(path)
#print(dir_list)

name = '%srlx'%(state_inp)
with open(name, 'w+') as inp:
    inp.truncate(0)

if 'GEOMETRY' not in dir_list:
    print('No GEOMETRY file is found')
    print('Written in %s'%(name))
    with open(name, 'a+') as inp:
        for line in nfinp:
            print('%s' %(line.rstrip("\n")),file = inp)    
else:
    print('Optimized cooridinate read from GEOMETRY')
    cps_opt   = np.zeros((katm,3),dtype = float)

    with open('GEOMETRY','r') as inp:
        geometry = inp.readlines()

    for i in range(katm):
        dummy = geometry[8+i]
        cps_opt[i,0:3] = dummy.split()[0:3]

    with open(name, 'a+') as inp:
        for i in range(start):
            print('%s' %(nfinp[i].rstrip("\n")),file = inp)   
        for i in range(katm):
            print('% 16.12f    % 16.12f    % 16.12f   %s %s %s'
                   %(cps_opt[i,0], cps_opt[i,1], cps_opt[i,2],nfinp[start+i].split()[3],nfinp[start+i].split()[4],nfinp[start+i].split()[5]),file = inp)   
        for i in range(start+katm,len(nfinp)):
            print('%s' %(nfinp[i].rstrip("\n")),file = inp)   
  
