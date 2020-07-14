from state import STATE_Input
from ase import Atom
from ase.io import read


#length is in angstrom 

# from ASE atoms
system = read('nfinp_1.xyz')
system.set_cell([(30,0,0),(0,30,0),(0,0,30)])

#from state
state = STATE_Input(system) 

#modify the input based on sub-class
state.kpoints.kpoints_mesh  = (1,1,1)
state.kpoints.kpoints_shift = (1,1,1)
state.pseudo.zeta = (0.1,0.0,0.1)
state.system.gmax = 4
state.system.gmaxp = 10
state.system.keg = 64


for i in range(3,5):
    print(i)
    state.constraints.imdtyp[i] = 0 
#state_test.kpoints_mesh  = (10,10,1)
#state_test.kpoints_shift = (2,2,1)

state.write_input ('nfinp_1')



#print(system.cell[:])

#print(atoms.get_chemical_formula())
#print(atoms.get_positions())
#print(atoms.get_positions())
#print(len(atoms.get_chemical_symbols()))
#print(len(atoms.get_chemical_symbols()))
#print(get_reduce_atom_list(atoms))
#print(get_ntyp(atoms))
#print(get_natm(atoms))
#a = cutoff() 
#print (gmax,gmaxp)
#print(get_ntyp())
#print(get_ntyp())
#len(atom)
#print(atom)i


