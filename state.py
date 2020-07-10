"""
Python3 script to generate the fixed input format of STATE code.
Atoms object from ASE is required.
Written by Thanh Ngoc Pham, OU.
"""

### import library
import numpy as np
### end 

### The unit is in default of fixed input format
ang2bohr = 1.8897259886
###

class STATE_Input:
    def __init__(self, atoms):
        """
        used several sub-classes
        """
        self.atoms       = atoms
        self.cell        = atoms.cell[:]     #print cell in 3 x 3 array in angstrom from ASE
        self.ntyp        = get_ntyp( atoms )     
        self.natm        = get_natm( atoms )
        self.space_group = space_group()
        self.control     = control()
        self.system      = system ()        #constains encut, smearing, spin, xctype, etc (try to duplicate the way QE works)
        self.kpoints     = kpoints()        #at present, only automatic kpoint is accepted
        self.electrons   = electrons ()    
        self.ions        = ions ()    
        self.pseudo      = pseudo (atoms)     
    def write_input (self,filename):
        write_state_input (self,filename)

class system:
    def __init__(self):
        self.gmax  = 6
        self.gmaxp = 20
        self.width = -0.002 
        self.xctype = 'ggapbe'
        self.kspin = 1
        self.nbztyp = 102
        self.ncord = 1
        self.ninv =  0

class control:
    def __init__(self):
        self.icond  = 0
        self.inipos = 0
        self.inivel = 0
        self.ininos = 0
        self.iniacc = 0
        self.cpumax = 86400 # one day
        self.ifstop =  0 

class electrons:
    def __init__(self):
        self.mixing_mode  = 6  # blugel
        self.mixing_what  = 1  # charge  
        self.iter_start   = 0
        self.nmd1        = 200 
        self.kbxmix       = 30
        self.mix_alpha    = 0.5
        self.nextst       = 1
        self.imsd         = 2  #   DAV is used as default, if rmm is used, nextst = 0 and imsd = 1

class ions:
    def __init__(self):
        self.nmd2  = 200  # blugel

class space_group:
    def __init__(self):
        self.num_space_group  = 1
        self.sg_type          = 0

class kpoints: 
    """  
    at present only automatic kpoint is accepted
    """
    kpoints_mesh  = np.zeros((1,3),dtype = int)
    kpoints_shift = np.zeros((1,3),dtype = int)


class pseudo:
    def __init__(self,atoms):
        self.iatomn = np.ones(get_ntyp(atoms),dtype = int)*29
        self.alfa  = np.ones(get_ntyp(atoms),dtype = float)* 0.5
        self.amion = np.ones(get_ntyp(atoms),dtype = float)*30  # mass of 30 amu is used as default
        self.iloc  = np.ones(get_ntyp(atoms),dtype = int)
        self.ivan  = np.ones(get_ntyp(atoms),dtype = int)
        self.zeta  = np.zeros(get_ntyp(atoms),dtype = float)

#def print_pseudo ( atoms, filename ):
#    atm_species =  get_reduce_atom_list(atoms)
#    """
#    print pseudo 
#    """
periodic = { 'H'   :'1'    ,
     'He'  :'2'    ,
     'B'   :'5'    ,
     'C'   :'6'    ,
     'N'   :'7'    ,
     'O'   :'8'    ,
     'F'   :'9'    ,
     'Ne'  :'10'   ,
     'Na'  :'11'   ,
     'Mg'  :'12'   ,
     'Al'  :'13'   ,
     'Si'  :'14'   ,
     'P'   :'15'   ,
     'S'   :'16'   ,
     'Cl'  :'17'   ,
     'Ar'  :'18'   ,
     'K'   :'19'   ,
     'Ca'  :'20'   ,
     'Sc'  :'21'   ,
     'Ti'  :'22'   ,
     'V'   :'23'   ,
     'Cr'  :'24'   ,
     'Mn'  :'25'   ,
     'Fe'  :'26'   ,
     'Co'  :'27'   ,
     'Ni'  :'28'   ,
     'Cu'  :'29'   ,
     'Zn'  :'30'   ,
     'Ga'  :'31'   ,
     'Ge'  :'32'   ,
     'As'  :'33'   ,
     'Se'  :'34'   ,
     'Br'  :'35'   ,
     'Kr'  :'36'   ,
     'Rb'  :'37'   ,
     'Sr'  :'38'   ,
     'Y'   :'39'   ,
     'Zr'  :'40'   ,
     'Nb'  :'41'   ,
     'Mo'  :'42'   ,
     'Tc'  :'43'   ,
     'Ru'  :'44'   ,
     'Rh'  :'45'   ,
     'Pd'  :'46'   ,
     'Ag'  :'47'   ,
     'Cd'  :'48'   ,
     'In'  :'49'   ,
     'Sn'  :'50'   ,
     'Ir'  :'77'   ,
     'Pt'  :'78'   ,
     'Au'  :'79'} 

def get_reduce_atom_list ( atoms ):
    """
    Get list of atomic symbol then reduce it.
    """
    ktyp_l =  atoms.get_chemical_symbols()
    if len ( ktyp_l ) > 1:
        for i in range (len( ktyp_l )-1, 0 , -1):
            for j in range ( i ):
                if ktyp_l [ j ] == ktyp_l [ i ]:
                    del ktyp_l [ i ]
                    break
    return ktyp_l

def get_ntyp ( atoms ):
    """
    Get number of atom type
    """
    list_atoms = get_reduce_atom_list( atoms )
    ntyp = len(list_atoms)
    return ntyp

def get_natm ( atoms ):
    """
    Get number of atom type
    """
    natm = len(atoms.get_positions())
    return natm

def print_cell ( cell,filename ):
    """
    print cell from numpy array
    """
    np.savetxt(filename,cell*ang2bohr, fmt = "% 16.12f") 

def print_cps ( atoms, filename ):
    """
    print atomic coordinate from numpy array in state format: cps_x cps_y cps_z iwei(=1) imdtyp ityp 
    """
    katm        =  get_natm(atoms)
    atm_species =  get_reduce_atom_list(atoms)
    cps         =  atoms.get_positions()
    ityp        =  np.ones( katm,dtype = int)
    imdtyp      =  np.ones( katm,dtype = int)
    ityp_list   =  atoms.get_chemical_symbols()

    for i in range(len(atm_species)):
        for j in range(katm):
            if ityp_list[j] == atm_species[i]: 
                ityp[j] = i+1

    for i in range(katm):
        print("% 16.12f   % 16.12f   % 16.12f %3s  %3s  %3s"
              %(cps[i,0],cps[i,1],cps[i,2],1,imdtyp[i], ityp[i])
              , file = filename)

def write_state_input (object, filename):
    with open(filename, 'w') as inp:
        print ('1 0 0 0 0 0',file = inp)
        print ( " % 4.2f % 4.2f %4i %4i %4i : gmax, gmaxp, ktyp, katm, katm2"
                 %(object.system.gmax, object.system.gmaxp,object.ntyp, object.natm, object.natm), file = inp)
        print ('%4i %4i :    num_space_group, type' 
                 %(object.space_group.num_space_group, object.space_group.sg_type),file = inp)
        print ('  Cartesian',file = inp)
        print_cell(object.cell,inp)
        print (' %4i %4i %4i %4i %4i %4i  : K_mesh  '
                 %(object.kpoints.kpoints_mesh[0],
                   object.kpoints.kpoints_mesh[1],
                   object.kpoints.kpoints_mesh[2],
                   object.kpoints.kpoints_shift[0],
                   object.kpoints.kpoints_shift[1],
                   object.kpoints.kpoints_shift[2]),file = inp)
        print ('%4i %4i :    num_space_group, type' 
                 %(object.system.ncord, object.system.ninv),file = inp)
        print_cps(object.atoms,inp)
        atm_species =  get_reduce_atom_list(object.atoms)
        for i in range(len(atm_species)):
            print('% 4.2f  % 4.2f  % 4.2f  % 2s   % 2s   % 4.2f '
                  %(int(periodic[atm_species[i]]), 
                   object.pseudo.alfa[i],
                   object.pseudo.amion[i],
                   object.pseudo.iloc[i],
                   object.pseudo.ivan[i],
                   object.pseudo.zeta[i] ),file=inp)
        print ('% 4i % 4i % 4i % 4i % 4i  : icond, inipos, inivel, ininos, iniacc'
                 %(object.control.icond,
                   object.control.inipos,
                   object.control.inivel,
                   object.control.ininos,
                   object.control.iniacc),file = inp)
        print ('   0    1   :  ipre ipri', file = inp)
        print ('%4i    %4i   0   %4i %4i :nmd1 nmd2 last_iter,cpumax,ifstop' 
                 %(object.electrons.nmd1,
                   object.control.cpumax,
                   object.ions.nmd2,
                   object.control.ifstop),file = inp)
