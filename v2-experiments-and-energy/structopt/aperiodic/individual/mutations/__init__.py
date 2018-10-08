import functools
import structopt.common.individual.mutations

from .move_atoms import move_atoms
from .move_surface_atoms import move_surface_atoms
from .move_surface_STEM import move_surface_STEM
from .move_atoms_group import move_atoms_group
from .rotate_cluster import rotate_cluster
from .twist import twist
from .swap_core_shell import swap_core_shell
from .rich2poor import rich2poor
from .poor2rich import poor2rich
from .flip_surface_atom import flip_surface_atom
from .permute_column_STEM import permute_column_STEM
from .permute_column_surface import permute_column_surface
from .permute_column_bulk import permute_column_bulk
from .move_surface_defects import move_surface_defects
from .enrich_surface import enrich_surface
from .enrich_bulk import enrich_bulk
from .enrich_surface_defects import enrich_surface_defects
from .enrich_surface_facets import enrich_surface_facets
from .permutation_STEM import permutation_STEM
from .increase_Z_STEM import increase_Z_STEM
from .decrease_Z_STEM import decrease_Z_STEM
from .enrich_surface_column import enrich_surface_column
from .enrich_bulk_column import enrich_bulk_column
from .rich2poor_column import rich2poor_column
from .poor2rich_column import poor2rich_column
from .move_column_defects import move_column_defects
from .move_column_random import move_column_random
from .add_atom_STEM import add_atom_STEM
from .add_atom_defects import add_atom_defects
from .add_atom_random import add_atom_random
from .remove_atom_STEM import remove_atom_STEM
from .remove_atom_defects import remove_atom_defects
from .remove_atom_random import remove_atom_random
from .move_surface_SCSA import move_surface_SCSA

move_surface_atoms.tag = 'MoSuAt'
move_surface_STEM.tag = 'MoSuSTEM'
move_atoms.tag = 'MoAt'
move_atoms_group.tag = 'MoAtGr'
rotate_cluster.tag = 'RoCl'
twist.tag = 'Twist'
swap_core_shell.tag = 'SwCoSh'
rich2poor.tag = 'Ri2Po'
poor2rich.tag = 'Pr2Ri'
permute_column_surface.tag = 'PeCoSu'
flip_surface_atom.tag = 'FlSuAt'
move_surface_defects.tag = 'MoSuDe'
enrich_surface.tag = 'EnSu'
enrich_bulk.tag = 'EnBu'
enrich_surface_defects.tag = 'EnSuDe'
enrich_surface_facets.tag = 'EnSuFa'
permutation_STEM.tag ='PeSTEM'
increase_Z_STEM.tag = 'InZSTEM'
decrease_Z_STEM.tag = 'DeZSTEM'
enrich_surface_column.tag = 'EnSuCo'
enrich_bulk_column.tag = 'EnBuCo'
rich2poor_column.tag = 'Ri2PoCo'
poor2rich_column.tag = 'Po2RiCo'
move_column_defects.tag = 'MoCoDe'
move_column_random.tag = 'MoCoRa'
add_atom_STEM.tag = 'AdAtSTEM'
add_atom_defects.tag = 'AdAtDe'
add_atom_random.tag = 'AdAtRa'
remove_atom_STEM.tag = 'ReAtSTEM'
remove_atom_defects.tag = 'ReAtDe'
remove_atom_random.tag = 'ReAtRa'
move_surface_SCSA.tag = 'MoSuSCSA'

class Mutations(structopt.common.individual.mutations.Mutations):

    @staticmethod
    @functools.wraps(move_surface_atoms)
    def move_surface_atoms(individual, max_natoms=0.2, move_CN=8, surf_CN=11):
        return move_surface_atoms(individual, max_natoms, move_CN, surf_CN)
    
    @staticmethod
    @functools.wraps(move_atoms)
    def move_atoms(individual, max_natoms=0.20):
        return move_atoms(individual, max_natoms)

    @staticmethod
    @functools.wraps(move_atoms_group)
    def move_atoms_group(individual, max_natoms=0.20):
        return move_atoms_group(individual, max_natoms)

    @staticmethod
    @functools.wraps(rotate_cluster)
    def rotate_cluster(individual, max_natoms=0.20):
        return rotate_cluster(individual, max_natoms)

    @staticmethod
    @functools.wraps(twist)
    def twist(individual, max_radius=0.90):
        return twist(individual, max_radius)
    
    @staticmethod
    @functools.wraps(swap_core_shell)
    def swap_core_shell(individual, max_natoms=0.2, surf_CN=11):
        return swap_core_shell(individual, max_natoms, surf_CN)

    @staticmethod
    @functools.wraps(move_surface_STEM)
    def move_surface_STEM(individual, STEM_parameters, move_CN=11, surf_CN=11,
                          filter_size=1, move_cutoff=0.5, surf_cutoff=0.5,
                          max_cutoff=0.5, min_cutoff=0.5):
        return move_surface_STEM(individual, STEM_parameters, move_CN, surf_CN,
                                 filter_size, move_cutoff, surf_cutoff,
                                 max_cutoff, min_cutoff)

    @staticmethod
    @functools.wraps(rich2poor)
    def rich2poor(individual):
        return rich2poor(individual)

    @staticmethod
    @functools.wraps(poor2rich)
    def poor2rich(individual):
        return poor2rich(individual)

    @staticmethod
    @functools.wraps(permute_column_surface)
    def permute_column_surface(individual, STEM_parameters, filter_size=0.5,
                               column_cutoff=0.5):
        return permute_column_surface(individual, STEM_parameters,
                                      filter_size, column_cutoff)

    @staticmethod
    @functools.wraps(flip_surface_atom)
    def flip_surface_atom(individual, surf_CN=11, cutoff=0.5):
        return flip_surface_atom(individual, surf_CN, cutoff)

    @staticmethod
    @functools.wraps(move_surface_defects)
    def move_surface_defects(individual, surf_CN=11):
        return move_surface_defects(individual, surf_CN)

    @staticmethod
    @functools.wraps(enrich_surface)
    def enrich_surface(individual, surf_CN=11, species=None):
        return enrich_surface(individual, surf_CN, species)

    @staticmethod
    @functools.wraps(enrich_bulk)
    def enrich_bulk(individual, surf_CN=11, species=None):
        return enrich_bulk(individual, surf_CN, species)

    @staticmethod
    @functools.wraps(enrich_surface_defects)
    def enrich_surface_defects(individual, surf_CN=11, species=None):
        return enrich_surface_defects(individual, surf_CN, species)

    @staticmethod
    @functools.wraps(enrich_surface_facets)
    def enrich_surface_facets(individual, surf_CN=11, species=None):
        return enrich_surface_facets(individual, surf_CN, species)

    @staticmethod
    @functools.wraps(permutation_STEM)
    def permutation_STEM(individual, STEM_parameters, filter_size=0.5,
                         move_cutoff=0.5, max_cutoff=0.5, min_cutoff=0.5):
        return permutation_STEM(individual, STEM_parameters, filter_size,
                                move_cutoff, max_cutoff, min_cutoff)

    @staticmethod
    @functools.wraps(increase_Z_STEM)
    def increase_Z_STEM(individual, STEM_parameters, filter_size=0.5,
                         move_cutoff=0.5, min_cutoff=0.5):
        return increase_Z_STEM(individual, STEM_parameters, filter_size,
                                move_cutoff, min_cutoff)

    @staticmethod
    @functools.wraps(decrease_Z_STEM)
    def decrease_Z_STEM(individual, STEM_parameters, filter_size=0.5,
                         move_cutoff=0.5, max_cutoff=0.5):
        return decrease_Z_STEM(individual, STEM_parameters, filter_size,
                                move_cutoff, max_cutoff)

    @staticmethod
    @functools.wraps(enrich_surface_column)
    def enrich_surface_column(individual, STEM_parameters, filter_size=0.5,
                              column_cutoff=0.5, species=None, surf_CN=11):
        return enrich_surface_column(individual, STEM_parameters, filter_size,
                                     column_cutoff, species, surf_CN)

    @staticmethod
    @functools.wraps(enrich_bulk_column)
    def enrich_bulk_column(individual, STEM_parameters, filter_size=0.5,
                           column_cutoff=0.5, species=None, surf_CN=11):
        return enrich_bulk_column(individual, STEM_parameters, filter_size,
                                  column_cutoff, species, surf_CN)
    
    @staticmethod
    @functools.wraps(rich2poor_column)
    def rich2poor_column(individual, STEM_parameters, filter_size=0.5,
                         column_cutoff=0.5, species=None, surf_CN=11):
        return rich2poor_column(individual, STEM_parameters, filter_size,
                                column_cutoff, species, surf_CN)

    @staticmethod
    @functools.wraps(poor2rich_column)
    def poor2rich_column(individual, STEM_parameters, filter_size=0.5,
                         column_cutoff=0.5, species=None, surf_CN=11):
        return poor2rich_column(individual, STEM_parameters, filter_size,
                                column_cutoff, species, surf_CN)

    @staticmethod
    @functools.wraps(move_column_defects)
    def move_column_defects(individual, cutoff=0.2, CN_factor=1.1):
        return move_column_defects(individual, cutoff, CN_factor)

    @staticmethod
    @functools.wraps(move_column_random)
    def move_column_random(individual, cutoff=0.2):
        return move_column_random(individual, cutoff)

    @staticmethod
    @functools.wraps(add_atom_STEM)
    def add_atom_STEM(individual, STEM_parameters, add_prob=None, permute=0.5,
                      filter_size=1, column_cutoff=0.2, surf_cutoff=0.5, min_cutoff=0.5):
        return add_atom_STEM(individual, STEM_parameters, add_prob, permute,
                             filter_size, column_cutoff, surf_cutoff, min_cutoff)

    @staticmethod
    @functools.wraps(add_atom_defects)
    def add_atom_defects(individual, add_prob=None, cutoff=0.2, CN_factor=1.1):
        return add_atom_defects(individual, add_prob, cutoff, CN_factor)

    @staticmethod
    @functools.wraps(add_atom_random)
    def add_atom_random(individual, add_prob=None, cutoff=0.2):
        return add_atom_random(individual, add_prob, cutoff)

    @staticmethod
    @functools.wraps(remove_atom_STEM)
    def remove_atom_STEM(individual, STEM_parameters, permute=True, remove_prob=None,
                         filter_size=1, remove_CN=11, remove_cutoff=0.5,
                         max_cutoff=0.5, column_cutoff=0.2):
        return remove_atom_STEM(individual, STEM_parameters, permute, remove_prob,
                                filter_size, remove_CN, remove_cutoff,
                                max_cutoff, column_cutoff)

    @staticmethod
    @functools.wraps(remove_atom_defects)
    def remove_atom_defects(individual, surf_CN=11):
        return remove_atom_defects(individual, surf_CN)

    @staticmethod
    @functools.wraps(remove_atom_random)
    def remove_atom_random(individual, surf_CN=11):
        return remove_atom_random(individual, surf_CN)
    
