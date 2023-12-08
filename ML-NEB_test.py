import os
from catlearn.optimize.mlneb import MLNEB
from ase.constraints import FixAtoms, FixedPlane
from ase.io import read
from ase.optimize import BFGS
from carmm.run.aims_path import set_aims_command
import numpy as np
from ase.neb import NEBTools
import matplotlib.pyplot as plt

set_aims_command(hpc="archer2", basis_set="light", defaults=2020)


def my_calc(periodic=True):
    # New method that gives a default calculator
    from carmm.run.aims_calculator import get_aims_and_sockets_calculator
    # On machines where ASE and FHI-aims are run separately (e.g. ASE on login node, FHI-aims on compute nodes)
    # we need to specifically state what the name of the login node is so the two packages can communicate
    sockets_calc, fhi_calc = get_aims_and_sockets_calculator(dimensions=2, k_grid=(5, 7, 1), #xc='pbe',
                                                             xc='libxc MGGA_X_MBEEF+GGA_C_PBE_SOL',
                                                             compute_forces="true")

    fhi_calc.set(xc_pre=['pbe', '10'],
                 # override_warning_libxc="true",
                 spin='none',
                 use_dipole_correction='True',
                 relativistic=('atomic_zora', 'scalar'),
                 compute_forces="true",
                 #many_body_dispersion_nl='beta=0.81',
                 charge_mix_param=0.02,
                 occupation_type="gaussian 0.05",
                 )

    # remove k_grid if not periodic
    if not periodic:
        fhi_calc.parameters.pop("k_grid")
        fhi_calc.parameters.pop("charge_mix_param")
        fhi_calc.parameters.pop("occupation_type")

    return sockets_calc


slab_initial = read("slab_initial.traj")
slab_final = read("slab_final.traj")

n = 7  # Desired number of images including start and end point, could be fraction
interpolation = 'idpp'

# Setup the Catlearn object for MLNEB
with my_calc() as calculator:
    neb_catlearn = MLNEB(start=slab_initial,
                         end=slab_final,
                         ase_calc=calculator,
                         n_images=n,
                         interpolation=interpolation,
                         restart=True)  # When True looks for evaluated_structures.traj

    # Run the NEB optimisation. Adjust fmax to desired convergence criteria, usually 0.01 ev/A
    neb_catlearn.run(fmax=0.01, trajectory='ML-NEB.traj', full_output=False, steps=75)

neb_tools = NEBTools(read("ML-NEB.traj@:"))

# Get the calculated barrier and the energy change of the reaction.
Ef, dE = neb_tools.get_barrier()

# Get the actual maximum force at this point in the simulation.
max_force = neb_tools.get_fmax()

# Create a figure like that coming from ASE-GUI.
fig = neb_tools.plot_band()
fig.savefig('barrier-mlneb.png')
print("Energy barrier:", Ef, " Delta E", dE, " max_force", max_force)

