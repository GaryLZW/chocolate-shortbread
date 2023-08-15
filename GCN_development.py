def coord_number(atoms, a=3.615, lattice='fcc'):
    # The list that stores coordination number for each atom
    cn_list = []
    # The list that stores the indices of first nearst neighbours for each atom
    # This will be a list of lists
    fnn_list = []

    if lattice == 'fcc':
        bond = round(a / 2 ** 0.5, 3)
    elif lattice == 'bcc':
        bond = round(a * (3 ** 0.5) / 2, 3)
    # Hi Luca, you can try to add an if statement to make this function also works for bcc

    # Distances with minimum image conversion.
    # A big enough model is still needed, e.g. (3*3*3)
    distances = atoms.get_all_distances(mic=True)

    for atom_i in atoms:
        i = atom_i.index
        cn = 0
        # List for first nearest neighbours of atom i
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

        # Append coordination number and first nearest neighbours to the lists every time the second j loop finishes.
        cn_list.append(cn)
        fnn_list.append(fnn)

    return cn_list, fnn_list


def general_coord_number(atoms, a, lattice, site):
    """
    :param atoms: Surface model (should be large enough, e.g.(3*3*3))
    :param a: lattice parameter
    :param lattice: crystal structure
    :param site: An atomic index for an ontop site or a list of atomic indices for a multi-atom adsorption site.
    :return: gcn. Generalized coordination number
    """
    # Hi Luca, you can try to complete this function for ontop sites using the coord_number function above.
    # The cn_max for fcc ontop sites is 12

    return gcn


from ase.build import bulk, fcc111

Cu = bulk('Cu', crystalstructure='fcc', a=3.615, cubic=True)
Cu = Cu.repeat((2, 2, 2))
Cu111 = fcc111('Cu', a=3.615, size=(3, 3, 3), vacuum=10, periodic=True)
gcn_fcc111 = general_coord_number(Cu111, a=3.615, lattice='fcc', site=[Cu111[-1].index])
print(gcn_fcc111)  # It should be 7.5
