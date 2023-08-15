from ase.build import fcc111, fcc110, fcc100, fcc211, add_adsorbate, molecule, surface, bulk
from ase import Atoms
from ase.constraints import FixAtoms, FixedLine
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.visualize import view

# Disabled. Won't work unless carmm is in your runtime environment.
# from software.run.aims_calculator import get_aims_calculator
# from software.run.aims_path import set_aims_command

emt = True

######ADSORBATE#########
## Define adsorbate, rotate and pre-optimize
molecule1 = molecule("CO")
# molecule1.rotate(0, 'x')

# molecule1 = Atoms('O')

if emt:
    molecule1.set_calculator(EMT())
# else:
#    set_aims_command()
#    molecule.set_calculator(get_aims_calculator("gas"))

## Optimize
molecule_opt = BFGS(molecule1)  # , trajectory="adsorbate.traj", restart="adsorbate.pckl")
molecule_opt.run(fmax=0.01)
e_opt_molecule = molecule1.get_potential_energy()

############SURFACE################################
from math import sqrt

## Edit these
# atomic_species='Cu'
# a_eos = 3.590
# unit_cell_depth=3
# unit_cell_width=3
# slab_depth=4
# vacuum_region_size=10.0

## Create surface
Cubulk = bulk('Cu', 'fcc', a=3.590, cubic=True)
slab = surface(Cubulk, (2, 1, 0), 4)
slab.center(vacuum=10.0, axis=2)
# slab = fcc111(atomic_species, a=a_eos,
# size=(unit_cell_width,unit_cell_depth,slab_depth),
# vacuum=vacuum_region_size)


## Set constraint for surface/bulk characteristics
#top_layers = 2  # leave x top layers relaxed
#mask0 = [atom.tag > top_layers for atom in slab]
#constraint0 = FixAtoms(mask=mask0)

#slab.set_constraint([constraint0])

if emt:
    slab.set_calculator(EMT())
# else:
#    slab.set_calculator(get_aims_calculator("periodic"))

## Preview your generated surface prior to running calculations
## Don't use this on HAWK or ISAMBARD! Job will end prematurely.
# from ase.visualize import view
# view(slab)
# view(molecule)

## Optimize the surface
surface_opt = BFGS(slab)  # , trajectory="surface.traj", restart='surface.pckl')
surface_opt.run(fmax=0.01)
e_opt_surface = slab.get_potential_energy()

## If happy with your structures, remove triple quotes to proceed

##########ADDING ADSORBATE ONTO THE SURFACE######
## Remove constraints and add adsorbate
slab.set_constraint()

## If build using ase.build can specify adsorption sites
add_adsorbate(slab, molecule1, 4.0, (slab[-1].x, slab[-1].y))
## Otherwise x-y coordinates - offset is specified in Angstrom
# add_adsorbate(slab, molecule, 2.0, position=(4.0, 2.4))

## Generate a new mask based on the changed number of atoms and constrain last two layers
#mask0 = [atom.tag > top_layers for atom in slab]
#constraint0 = FixAtoms(mask=mask0)

# Fix all the Cu atoms
mask1 = [atom.symbol == 'Cu' for atom in slab]
constraint0 = FixAtoms(mask=mask1)
# Uncomment these lines to constrain the relaxation of adsorbate to the z axis
indices1 = [atom.index for atom in slab if atom.symbol == 'C' or atom.symbol == 'O']
constraint1 = FixedLine(indices1[0], direction=[0, 0, 1])
if len(molecule1) == 2:
    constraint2 = FixedLine(indices1[1], direction=[0, 0, 1])  # if using oxygen, commennt constraint2
    slab.set_constraint([constraint0, constraint1, constraint2])
else:
    slab.set_constraint([constraint0, constraint1])

# Add constaint1 to set constraint, i.e. change [constraint0] to [constraint0, constraint1]


if emt:
    slab.set_calculator(EMT())
# else:
#    slab.set_calculator(get_aims_calculator("periodic"))

## Optimize
dyn = BFGS(slab)  # , trajectory='ads_slab.traj', restart="ads_slab.pckl")
dyn.run(fmax=0.01)  # tighten to min 0.01 eV/A for actual calculation
e_opt_slab = slab.get_potential_energy()

print("Energy of adsorbate: ", e_opt_molecule)
print("Energy of a clean surface: ", e_opt_surface)
print("Energy of an optimised slab: ", e_opt_slab)
Eb = (e_opt_slab - (e_opt_molecule + e_opt_surface))
print("Binding Energy: ", Eb)

# assert(abs(Eb - -0.173372) < 1e-5)

view(slab)
