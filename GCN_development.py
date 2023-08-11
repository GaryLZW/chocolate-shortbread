def coord_number(atoms, a=3.615, lattice='fcc'):
   from ase import Atoms
   cn = []
   if lattice == 'fcc':
       bond = a / 2 ** 0.5
   distances = atoms.get_all_distances()
   CN = 0
   for atom_i in atoms:
       i = atom_i.index
       for atom_j in atoms:
           j = atom_j.index
           if i == j:
               continue
           if distances[i][j] == bond:
               CN += 1
   cn.append(bond)
   return CN
