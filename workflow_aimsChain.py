import aimsChain
import ase.calculators.aims
import scipy
import numpy as np
from ase.io.aims import write_aims

from carmm.run.aims_path import set_aims_command
from carmm.run.aims_calculator import get_aims_calculator
from ase.io import read
import os



def write_chain_control(parameters: dict):
    lim = '#' * 80
    chain_in = open(file='chain.in', mode='w')
    keyword_category = {'Environment setting': ['run_aims', 'initial_file', 'final_file'],
                        'General settings for task': ['method', 'n_images', 'force_thres', 'optimizer', 'restart',
                                                      'neb_spring_constant', 'global_optimizer'],
                        'External starting geometries': ['external_geometry', 'resample'],
                        'Climbing image setting': ['use_climb', 'climb_thres', 'interpolated_climb', 'climb_mode',
                                                   'climb_global_optimizer', 'climb_optimizer', 'climb_control',
                                                   'lbfgs_alpha', 'lbfgs_memory', 'lbfgs_maxstep', 'bfgs_alpha',
                                                   'bfgs_maxstep', 'fire_dt'],
                        'Growing string method': ['use_gs_method', 'gs_n_images', 'gs_thres', 'gs_optimizer',
                                                  'gs_global_optimizer', 'fire_maxstep'],
                        'Periodic system setting': ['periodic_interpolation', 'xyz_lattice', 'map_unit_cell']}
    for category, keywords in keyword_category.items():
        chain_in.write(lim + '\n' + '##' + category + '\n' + lim + '\n')
        [chain_in.write('%-20s%s\n' % (key, value)) for key, value in parameters.items()
         if key in keywords]
    return 0


class AimsChainWorkFlow:
    def __init__(self, fhi_calc: ase.calculators.aims.Aims, initial_struc: ase.Atoms, final_struc: ase.Atoms,
                 geo_constrain: bool = False, **kwargs):
        self.parameters = dict(kwargs)
        self.initial = initial_struc
        self.final = final_struc
        self.calc = fhi_calc
        self.geo_constrain = geo_constrain

        self.parameters['run_aims'] = self.calc.command

    def write_chain_inputs(self):
        write_chain_control(self.parameters)
        write_aims(os.path.join(self.calc.directory, self.parameters['initial_file']), atoms=self.initial,
                   geo_constrain=self.geo_constrain)
        write_aims(os.path.join(self.calc.directory, self.parameters['final_file']), atoms=self.final,
                   geo_constrain=self.geo_constrain)
        self.calc.write_control(atoms=self.initial, filename='control.in')


set_aims_command(hpc='hawk', basis_set='light', defaults=2020)
calc = get_aims_calculator(dimensions=2, k_grid=(5, 7, 1), xc='libxc MGGA_X_MBEEF+GGA_C_PBE_SOL')
calc.set(xc_pre=['pbe', '10'],
         spin='none',
         use_dipole_correction='True',
         relativistic=('atomic_zora', 'scalar'),
         compute_forces="true",
         charge_mix_param=0.02,
         occupation_type="gaussian 0.05",
         )
initial = read("initial.traj")
final = read("final.traj")
string_chain = AimsChainWorkFlow(fhi_calc=calc, initial_struc=initial, final_struc=final, initial_file='ini.in',
                                 final_file='fin.in', method='string', n_images=7, force_thres=0.01, optimizer='BFGS',
                                 use_climb='true', climb_thres=0.05)
print(string_chain.parameters['initial_file'])
string_chain.write_chain_inputs()

