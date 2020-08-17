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


def read_statedos(nfout,nfinp):
    state_out = nfout
    state_in  = nfinp
    AO_LDOS   =  'AO_LDOS:'
    kspin     =  'kspin'
    npdos     =  'npdosao'
    index     = {AO_LDOS :[] }
    index_in  = {kspin   :[],
                 npdos   :[] }

###Read AOLDOS from nfout
    with open(state_out, 'r') as nfout:
        nfout_lines = nfout.readlines()
        for idx,line in enumerate(nfout_lines):
            for identifier in index:
                if identifier in line:
                    index[identifier].append(idx)


###Read input for AO_LDOS (kspin, nodosao and atom_index in AOLDOS)
    with open(state_in, 'r') as nfinp:
        nfinp_lines = nfinp.readlines()
        for idx,line in enumerate(nfinp_lines):
            for identifier in index_in:
                if identifier in line:
                    index_in[identifier].append(idx)

    kspin  = int(nfinp_lines[index_in[kspin][-1]].split()[1])
    if kspin == 1:
        print('Non spin polarization calculation is performed')
    else:
        print('Spin polarization calculation is performed')

    npdosao = int(nfinp_lines[index_in[npdos][-1]].split()[0])
    print(int(nfinp_lines[index_in[npdos][-1]].split()[0]))
    print('Number of atoms in AO_LDOS is %4d' %(npdosao))

    atom_index = np.zeros(npdosao,dtype = int)

    for i in range(len(atom_index)):
        atom_index[i] = int(nfinp_lines[index_in[npdos][-1]+i+1].split()[0])
    print('List of atom index in AO_LDOS:')
    print(atom_index)

    natm = int(nfout_lines[index[AO_LDOS][-1]].split()[1])
    npdose = int((len(index[AO_LDOS])-1)/npdosao/kspin) 


### Write the output
### For non spin calculation
    if kspin == 1:
        for i in atom_index:
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
## e, e,s,p,d,tdos,px,py,pz,dz2,dxxyy,dxy,dyz,dzx
            with open(filename+'.npy','wb') as f:
                np.save(f, dname)
    else:
###For spin polarization calculation
        for i in atom_index:
            spinup   = 'pdos_up_%04d.dat' %(i)
            spindown = 'pdos_down_%04d.dat' %(i)
            with open(spinup, 'w') as dos:
                for j in range(1,(len(index[AO_LDOS]))):
                    atm_index = int((nfout_lines[index[AO_LDOS][j]].split()[1]))
                    spin      = int((nfout_lines[index[AO_LDOS][j]].split()[2]))
                    if (atm_index == i) & (spin ==2):
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
        for i in atom_index:
            spinup = 'pdos_up_%04d' %(i)
            dname  = spinup
            dname  = np.loadtxt(fname=spinup+'.dat',
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
## e, e,s,p,d,tdos,px,py,pz,dz2,dxxyy,dxy,dyz,dzx
            with open(spinup+'.npy','wb') as f:
                np.save(f, dname)
# for down spin 
            with open(spindown, 'w') as dos:
                for j in range(1,(len(index[AO_LDOS]))):
                    atm_index = int((nfout_lines[index[AO_LDOS][j]].split()[1]))
                    spin      = int((nfout_lines[index[AO_LDOS][j]].split()[2]))
                    if (atm_index == i) & (spin ==1):
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
        for i in atom_index:
            spindown = 'pdos_down_%04d' %(i)
            dname  = spindown
            dname  = np.loadtxt(fname=spindown+'.dat',
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
## e, e,s,p,d,tdos,px,py,pz,dz2,dxxyy,dxy,dyz,dzx
            with open(spindown+'.npy','wb') as f:
                np.save(f, dname)
