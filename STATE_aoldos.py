"""
Python3 script to read atomic orbital projected local density of state (AO_LDOS) result from STATE code

AO_LDOS written in 'pdos_%04d.dat' %(i) file (i is atom index) has following format:
 e,s,p,d,tdos,px,py,pz,dz2,dxxyy,dxy,dyz,dzx

To load the result and plot using mathplotlib in python3,
    First, load the .npy
        atom1 = np.load('pdos_0001.npy')
    Secondly, plot by call the component 
    plt.plot(atom2['e'], atom2['tot'], '-', markersize = 10,linewidth=3).
    The component name is similar to above written for .dat files

Written by Thanh Ngoc Pham, OU.
"""

### import library
import numpy as np
### end 

### The unit is in default of fixed input format
ang2bohr = 1.8897259886
###


def read_statedos(nfout):
    filename = nfout
    AO_LDOS   =  'AO_LDOS:'
    index = {AO_LDOS :[] }

    with open(filename, 'r') as nfout:
        nfout_lines = nfout.readlines()
        for idx,line in enumerate(nfout_lines):
            for identifier in index:
                if identifier in line:
                    index[identifier].append(idx)

    natm = int(nfout_lines[index[AO_LDOS][-1]].split()[1])
    npdose = int((len(index[AO_LDOS])-1)/natm) 

    ### write the dat files of AO_LDOS of each atom

    for i in range(1,natm+1):
        filename = 'pdos_%04d.dat' %(i)
        with open(filename, 'w') as dos:
            for j in range(1,(len(index[AO_LDOS]))):
                atm_index = int((nfout_lines[index[AO_LDOS][j]].split()[1]))
                if atm_index == i:
                  energy = float(nfout_lines[index[AO_LDOS][j]].split()[3])
                  s      = float(nfout_lines[index[AO_LDOS][j]].split()[4])
                  px     = float(nfout_lines[index[AO_LDOS][j]].split()[5])
                  py     = float(nfout_lines[index[AO_LDOS][j]].split()[6])
                  pz     = float(nfout_lines[index[AO_LDOS][j]].split()[7])
                  dz2    = float(nfout_lines[index[AO_LDOS][j]].split()[8])
                  dxxyy  = float(nfout_lines[index[AO_LDOS][j]].split()[9])
                  dxy    = float(nfout_lines[index[AO_LDOS][j]].split()[10])
                  dyz    = float(nfout_lines[index[AO_LDOS][j]].split()[11])
                  dzx    = float(nfout_lines[index[AO_LDOS][j]].split()[12])
                  p = px+py+pz
                  d = dz2+dxxyy+dxy+dyz+dzx
                  totdos = float(nfout_lines[index[AO_LDOS][j]].split()[13])
                  print('%+8.5f    %+8.5f    %+8.5f   %+8.5f    %+8.5f    %+8.5f    %+8.5f   %+8.5f  %+8.5f    %+8.5f    %+8.5f    %+8.5f    %+8.5f'
                         %(energy,  s,       p,        d,        totdos,   px,       py,      pz,     dz2,      dxxyy,      dxy,      dyz,      dzx),
                         file = dos)



    ### store numpy binary file
    for i in range(1,natm+1):
        filename = 'pdos_%04d' %(i)
        dname  = filename
        dname  = np.loadtxt(fname=filename+'.dat',
                     dtype =([('e'     , float),
                             ( 's'     , float),
                             ( 'p'     , float),
                             ( 'd'     , float),
                             ( 'tdos'  , float),
                             ( 'px'    , float),
                             ( 'py'    , float),
                             ( 'pz'    , float),
                             ( 'dz2'   , float),
                             ( 'dxxyy' , float),
                             ( 'dxy'   , float),
                             ( 'dyz'   , float),
                             ( 'dzx'   , float)]))

# e,s,p,d,tdos,px,py,pz,dz2,dxxyy,dxy,dyz,dzx

     
        with open(filename+'.npy','wb') as f:
            np.save(f, dname)
