import aimsChain
import ase.calculators.aims
import scipy
import numpy as np
from ase.io.aims import write_aims


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
        [chain_in.write('%-35s%s\n' % (key, value)) for key, value in parameters.items()
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

    def write_chain_inputs(self):
        write_chain_control(self.parameters)
        write_aims(fd=self.parameters['initial_file'], atoms=self.initial, geo_constrain=self.geo_constrain)
        write_aims(fd=self.parameters['final_file'], atoms=self.final, geo_constrain=self.geo_constrain)
        self.calc.write_control(atoms=self.initial, filename='control.in')
