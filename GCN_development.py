def coord_number(atoms, a=3.615, lattice='fcc'):
   #List that stores coordination number for each atom
   cn_list = []
   
   #Possible error if have different decimal points?
   #Should we use the smallest non-zero value from distances?
   if lattice == 'fcc':
       bond = a / 2 ** 0.5
   distances = atoms.get_all_distances()
   
   
   #CN = 0 This line should be at the start of the first i loop as
   #we should count coordination number from zero for each atom
   for atom_i in atoms:
      i = atom_i.index
      cn = 0
      #Counting coordination number for atom i
      for atom_j in atoms:
         j = atom_j.index
         #Skip the iteration if we are considering the distance 
         #between atom i and itself
         if i == j:
            continue
         #Check if atom i and atom j are first nearst neighbours
         if distances[i][j] == bond:
            cn += 1
      #Append coordination number to the list every time the second j loop finishes.  
      cn_list.append(cn)
   
   return cn_list
