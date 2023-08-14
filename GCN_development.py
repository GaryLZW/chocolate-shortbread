def coord_number(atoms, a=3.615, lattice='fcc'):
    # List that stores coordination number for each atom
    cn_list = []
    # List that stores the indices of first nearst neighbours for each atom
    fnn_list = []
    # Possible error if have different decimal points?
    # Should we use the smallest non-zero value from distances?
    if lattice == 'fcc':
        bond = round(a / 2 ** 0.5, 3)

    # Distances with minimum image conversion.
    # A big enough model is still needed.
    distances = atoms.get_all_distances(mic=True)

    for atom_i in atoms:
        i = atom_i.index
        cn = 0
        fnn = []
        # Counting coordination number for atom i
        for atom_j in atoms:
            j = atom_j.index
            # Skip the iteration if we are considering the distance
            # between atom i and itself
            if i == j:
                continue
            # Check if atom i and atom j are first nearst neighbours
            if round(distances[i][j], 3) == bond:
                cn += 1
                fnn.append(j)
        # Append coordination number to the list every time the second j loop finishes.
        cn_list.append(cn)
        fnn_list.append(fnn)

    return cn_list, fnn_list


from ase.build import bulk, fcc111

Cu = bulk('Cu', crystalstructure='fcc', a=3.615, cubic=True)
Cu = Cu.repeat((2, 2, 2))
Cu111 = fcc111('Cu', a=3.615, size=(3, 3, 3), vacuum=10, periodic=True)
cn_list, fnn_list = coord_number(Cu111, a=3.615, lattice='fcc')
print(cn_list)
