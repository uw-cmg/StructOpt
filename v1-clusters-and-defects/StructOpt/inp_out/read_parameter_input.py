import os
import copy
import random
import time
import numpy
import logging
import json
from StructOpt.inp_out.loggerUtils import initialize_logger
from StructOpt.tools.check_atomlist_concentration import check_atomlist_concentration
import pdb
try:
    from mpi4py import MPI
except:
    pass

def read_parameter_input(input, logger):
    """Function to convert input string, file, or dictionary to a dictionary to contain 
        the parameters for use by the optimizer class.
        input:
            input : Can be dictionary, string, or filename
            logfile : Name of logfile to write to. Default is None
        output:
            parameters : dictionary for defining parameters for optimizer with defaults
    """
    parameters = None
    if isinstance(input,dict):
        # If supplied input is already a dictionary set parameters equal to that dictionary
        parameters=input
    else:
        #Check to see if input is a filename
        if os.path.exists(input):
            parameters = json.load(open(input))
            # parameters = eval("%s"%open(input).read()) # issues with True/true
    if 'parallel' not in parameters:
        parameters['parallel'] = False
    if 'filename' not in parameters:
        parameters['filename']='Output'
    try:
        rank = MPI.COMM_WORLD.Get_rank()
    except:
        rank = 0
    if logger:
        if 'loggername' not in parameters:
            parameters['loggername']= '{0}-rank{1}-{2}.log'.format(parameters['filename'],rank,time.strftime("%Y_%m%d_%H%M%S"))
            logger = initialize_logger(parameters['loggername'])
    else:
        parameters['loggername'] = None
        logger = dummy_logger_no_write()
    
    if 'modules' not in parameters:
        logger.warning('modules not set.  Default values set to LAMMPS.\n')
        parameters['modules'] = ['LAMMPS']
    if 'relaxation' not in parameters:
        if 'LAMMPS' in parameters['modules']:
            parameters['relaxation'] = 'LAMMPS'
        elif 'VASP' in parameters['modules']:
            parameters['relaxation'] = 'VASP'
        else:
            parameters['relaxation'] = None
    if 'weights' not in parameters:
        parameters['weights'] = [1.0 for m in range(len(parameters['modules']))]

    if ('atomlist' not in parameters):
        #Stop program if atom list parameter not in input
        logger.critical("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx',Concentration,Mass,Chemical Potential)]")
        logger.critical("Current parameters include:\n" + repr(parameters))
        raise RuntimeError("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx',Concentration,Mass,Chemical Potential)]")
    else:
        if isinstance(parameters['atomlist'],str):
            parameters['atomlist'] = eval(parameters['atomlist'].strip())
        if not isinstance(parameters['atomlist'],list):
            logger.critical('Something is wrong with atomlist parameter: {0}'.format(parameters['atomlist']))
            raise RuntimeError("Input file/string/dictionary must include an atomlist defined as 'atomlist':[('Xx',Concentration,Mass,Chemical Potential)]: {0}".format(parameters['atomlist']))
    if 'structure' not in parameters:
        #Stop program if structure parameter not in input
        logger.critical("Input file/dictionary must include a structure for the simulation as 'structure':'Cluster/Crystal/Defect'")
        logger.debug("Current parameters include:\n"+repr(parameters))
        raise RuntimeError("Input file/dictionary must include a structure for the simulation as 'structure':'Cluster/Crystal/Defect'")
    for one in parameters['atomlist']:
        if len(one) != 4:
            #Stop program if atom list parameter not properly formatted
            logger.critical('Format of atom list not correct. Must be [(symbol,concentration,mass,potential)]')
            logger.debug('Issue in section : {0}'.format(one))
            logger.debug('Current atomlist is formatted as : {0}'.format(parameters['atomlist']))
            raise RuntimeError('Format of atom list not correct. Must be [(symbol,concentration,mass,potential)]')
    if 'natoms' not in parameters:
        parameters['natoms'] = int(sum([abs(c) for ind,c,m,u in parameters['atomlist']]))
        if rank==0:
            logger.warning('Number of atoms in simulation not set')
            logger.warning('Assuming natoms = {0}'.format(parameters['natoms']))
    parameters['atomlist'] = check_atomlist_concentration(parameters['atomlist'],parameters['natoms'],parameters['loggername'])
    if 'optimizer_type' not in parameters:
        parameters['optimizer_type'] = 'GA'
        if rank==0:
            logger.info('optimizer_type not set.  Default values set to GA.')
    if parameters['optimizer_type'] == 'GA':
        nindiv = 10
        genealogy = True
        nbests = 100
        algtype = 'lambda+mu'
        cxpb = 0.8
        mutpb = 0.15
        natselectscheme = 'tournament'
    else:
        nindiv = 1
        genealogy = False
        nbests = 100
        algtype = parameters['optimizer_type']
        cxpb = 0.0
        mutpb = 1.0
    if parameters['optimizer_type'] == 'SA':
        natselectscheme = 'metropolis'
        predator = 'adapting'
    elif parameters['optimizer_type'] == 'BH':
        natselectscheme = 'metropolis'
    else:
        natselectscheme = 'best'
    if 'nindiv' not in parameters:
        parameters['nindiv'] = nindiv
        if rank==0:
            logger.info('Setting number of individuals in population (nindiv) = {0}'.format(parameters['nindiv'])) 
    
    #Parameters for output
    if 'genealogy' not in parameters:
        parameters['genealogy'] = genealogy
        if rank==0:
            logger.info('Setting genealogy = {0}'.format(parameters['genealogy']))
    if 'output_format' not in parameters:
        parameters['output_format'] = 'fitness'
        if rank==0:
            logger.info('Setting output format = {0}'.format(parameters['output_format']))
    if 'allenergyfile' not in parameters:
        parameters['allenergyfile'] = False
        if rank==0:
            logger.info('Setting allenergyfile = {0}'.format(parameters['allenergyfile']))
    if 'best_inds_list' not in parameters:
        parameters['best_inds_list'] = True
        if rank==0:
            logger.info('Setting best_inds_list = {0}'.format(parameters['best_inds_list']))
    if 'number_of_bests' not in parameters:
        parameters['number_of_bests'] = nbests
        if rank==0:
            logger.info('Setting number_of_bests = {0}'.format(parameters['number_of_bests']))
    if 'indiv_defect_write' not in parameters:
        parameters['indiv_defect_write'] = False
        if rank==0:
            logger.info('Setting indiv_defect_write = {0}'.format(parameters['indiv_defect_write']))
    if 'vacancy_output' not in parameters:
        parameters['vacancy_output'] = False
        if rank==0:
            logger.info('Setting vacancy_output = {0}'.format(parameters['vacancy_output']))
    if 'restart_optimizer' not in parameters:
        parameters['restart_optimizer'] = False
    #Parameters for post-processing
    if 'lattice_concentration' not in parameters:
        parameters['lattice_concentration'] = False
        if rank==0:
            logger.info('Setting lattice_concentration = {0}'.format(parameters['lattice_concentration']))
    if 'postprocessing' not in parameters:
        parameters['postprocessing'] = False
        if rank==0:
            logger.info('Setting postprocessing = {0}'.format(parameters['postprocessing']))
    if 'genealogytree' not in parameters:
        parameters['genealogytree'] = False
        if rank==0:
            logger.info('Setting genealogytree = {0}'.format(parameters['genealogytree']))
    
    #Parameters for general algorithm
    if 'seed' not in parameters:
        parameters['seed']=random.randint(0,10)
        if rank==0:
            logger.info('Setting Random number seed (seed) to {0}'.format(parameters['seed']))
    if 'forcing' not in parameters:
        parameters['forcing'] = 'Concentration'
        if rank==0:
            logger.info('Setting forcing = {0}'.format(parameters['forcing']))
            logger.info('Assuming forcing concentration control')
    if 'debug' not in parameters:
        parameters['debug'] = ['None']
        if rank==0:
            logger.info('Setting debug = {0}'.format(parameters['debug']))
        if 'None' not in parameters['debug']: 
            print '***** DEBUGGING RUN *****'
    if 'algorithm_type' not in parameters:
        parameters['algorithm_type'] = algtype
        if rank==0:
            logger.info('Setting algorithm type = {0}'.format(parameters['algorithm_type']))
    if 'migration_intervals' not in parameters:
        parameters['migration_intervals'] = 5
        if rank==0:
            logger.info('Setting migration_intervals = '.format(parameters['migration_intervals']))
    if 'migration_percent' not in parameters:
        parameters['migration_percent'] = 0.05
        if rank==0:
            logger.info('Setting migration_percent = {0}'.format(parameters['migration_percent']))
    if 'fingerprinting' not in parameters:
        parameters['fingerprinting'] = False
        if rank==0:
            logger.info('Setting fingerprinting = {0}'.format(parameters['fingerprinting']))
    if 'fpbin' not in parameters:
        parameters['fpbin'] = 0.25
        if rank==0:
            logger.info('Setting fingerprint bin to {0}'.format(parameters['fpbin']))
    if 'fpcutoff' not in parameters:
        parameters['fpcutoff'] = 15.0
        if rank==0:
            logger.info('Setting fingerprint cutoff distance to {0}'.format(parameters['fpcutoff']))
    parameters['bulkfp'] = None
    if 'fixed_region' not in parameters:
        parameters['fixed_region'] = False
        if rank==0:
            logger.info('Setting fixed_region = {0}'.format(parameters['fixed_region']))
    if 'rattle_atoms' not in parameters:
        parameters['rattle_atoms'] = False
        if rank==0:
            logger.info('Setting rattle_atoms = {0}'.format(parameters['rattle_atoms']))
    if 'constrain_position' not in parameters:
        parameters['constrain_position'] = False
        if rank==0:
            logger.info('Setting constrain_position = {0}'.format(parameters['constrain_position']))
    if 'restart' not in parameters:
        parameters['restart'] = False
        if rank==0:
            logger.info('Setting restart = {0}'.format(parameters['restart']))
    if 'restart_ints' not in parameters:
        parameters['restart_ints'] = 0
        if rank==0:
            if parameters['restart']:
                logger.info('Setting restart_ints = {0}'.format(parameters['restart_ints']))
    
    # Parameters to generate the population and individual
    if 'r_ab' not in parameters:
        parameters['r_ab'] = 2.5
        if rank==0:
            logger.info('Setting r_ab = {0}'.format(parameters['r_ab']))
    if 'size' not in parameters:
        parameters['size'] = parameters['natoms']**0.33333*parameters['r_ab']
        if rank==0:
            logger.info('Setting size to (natoms)^(1/3)*r_ab = {0}'.format(parameters['size']))
    if 'generate_flag' not in parameters:
        parameters['generate_flag']='box'
        if rank==0:
            logger.info('Setting default generation scheme = {0}'.format(parameters['generate_flag']))
    parameters['solidbulk'] = None
    if 'sf' not in parameters:
        parameters['sf'] = 1.75
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting size factor for Defect (sf) = {0}'.format(parameters['sf']))
    if 'supercell' not in parameters:
        parameters['supercell'] = (1,1,1)
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting supercell for Defect (supercell) = {0}'.format(parameters['supercell']))
    if 'solidfile' not in parameters:
        if parameters['structure'] == 'Defect':
            logger.critical('Must provide a file for bulk solid if running a defect simulation.')
            raise RuntimeError('Error: Bulk for Defect not specified. Enter name of file for bulk structure as solidfile parameter')
        else:
            parameters['solidfile'] = None
    if 'solidcell' not in parameters:
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.warning('Warning cell size for Bulk Solid not specified assuming distance between 1st and 3rd atom')
        parameters['solidcell'] = None
    if 'evalsolid' not in parameters:
        parameters['evalsolid'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Not evaluating Solid')
    if 'finddefects' not in parameters:
        parameters['finddefects'] = True
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting finddefects = {0}'.format(parameters['finddefects']))
    if 'trackvacs' not in parameters:
        parameters['trackvacs'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting trackvacs = {0}'.format(parameters['trackvacs']))
    if 'trackswaps' not in parameters:
        parameters['trackswaps'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting trackswaps = {0}'.format(parameters['trackswaps']))
    if 'random_loc_start' not in parameters:
        parameters['random_loc_start'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting random_loc_start = {0}'.format(parameters['random_loc_start']))
    if 'random_vac_start' not in parameters:
        parameters['random_vac_start'] = False
        if rank==0:
            if parameters['structure'] == 'Defect':
                logger.info('Setting random_vac_start = {0}'.format(parameters['random_vac_start']))
    if 'purebulkenpa' not in parameters:
        if parameters['structure'] == 'Defect':
            parameters['purebulkenpa'] = None
    if 'natomsbulk' not in parameters:
        if parameters['structure'] == 'Defect':
            parameters['natomsbulk'] = None
    if 'surfacefile' not in parameters:
        if parameters['structure'] == 'Surface':
            logger.critical('Must provide a file for bulk solid if running a defect simulation.')
            raise RuntimeError('Error: Bulk for Defect not specified. Enter name of file for bulk structure as solidfile parameter')
        else:
            parameters['surfacefile'] = None
    if 'surfacecell' not in parameters:
        parameters['surfacecell'] = None
    if 'surftopthick' not in parameters:
        parameters['surftopthick'] = 0        
    if 'cell_shape_options' not in parameters:
        parameters['cell_shape_options'] = ['cubic', 'hexagonal', 'triclinic', \
        'monoclinic', 'orthorhombic', 'tetragonal']
        if rank==0:
            if parameters['structure'] == 'Crystal':
                logger.info('Assuming following cell shape options: {0}'.format(parameters['cell_shape_options']))
    if 'alloy' not in parameters:
        parameters['alloy'] = True
        if rank==0:
            logger.info('Setting alloy = {0}'.format(parameters['alloy']))
    
    if 'large_box_size' not in parameters:
        parameters['large_box_size']=500.0
        if rank == 0:
            if parameters['structure']=='Cluster':
                logger.info('Setting large_box_size to {0}'.format(parameters['large_box_size']))
    
    #Parameters for Crossovers
    if 'cxpb' not in parameters:
        parameters['cxpb'] = cxpb
        if rank == 0:
            logger.info('Setting crossover probability (cxpb) = {0}'.format(parameters['cxpb']))
    if 'cx_scheme' not in parameters:
        parameters['cx_scheme'] = 'cxtp'
        if rank == 0:
            logger.info('Assuming two-point crossover.  Setting cx_scheme = {0}'.format(parameters['cx_scheme']))
    if 'selection_scheme' not in parameters:
        parameters['selection_scheme'] = 'tournament2'
        if rank == 0:
            logger.info('Setting selection_scheme = {0}'.format(parameters['selection_scheme']))
    
    #Parameters for Mutations
    if 'mutpb' not in parameters:
        parameters['mutpb'] = mutpb
        if rank == 0:
            logger.info('Setting mutation probability (mutpb) = {0}'.format(parameters['mutpb']))
    if 'mutation_options' not in parameters:
        if parameters['structure']=='Cluster':
            parameters['mutation_options']=['lattice_alteration','rotation',\
            'permutation','scale_size']
        elif parameters['structure']=='Crystal':
            parameters['mutation_options']=['lattice_alteration','rotation',\
            'permutation','scale_size', 'cell_shape', 'lammps_box_relax']
        elif parameters['structure']=='Defect':
            parameters['mutation_options']=['lattice_alteration','rotation','permutation']
        if rank == 0:
            logger.info('Setting mutations options = {0}'.format(parameters['mutation_options']))
    BHFlag=False
    for one in parameters['mutation_options']:
        if 'basin_hop' in one:
            BHFlag=True
    if 'bh_steps' not in parameters:
        parameters['bh_steps']=100
        if rank == 0:
            if BHFlag==True:
                logger.warning('Max steps not specified for Basin Hop mutation, setting bh_steps = {0}'.format(parameters['bh_steps']))
    if 'bh_temp' not in parameters:
        parameters['bh_temp'] = 1000*8.617385692256675e-05
        if rank == 0:
            if BHFlag==True:
                logger.warning('Temperature not set for Basin Hop mutation, setting bh_temp in kT = {0}'.format(parameters['bh_temp']))
    if 'mutant_add' not in parameters:
        parameters['mutant_add'] = False
        if rank == 0:
            logger.info('Setting mutant_add = {0}'.format(parameters['mutant_add']))
    if 'quench_max_temp' not in parameters:
        parameters['quench_max_temp'] = 1000
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Setting quench_max_temp = {0}'.format(parameters['quench_max_temp']))
    if 'quench_min_temp' not in parameters:
        parameters['quench_min_temp'] = 2
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Seting quench_min_temp = {0}'.format(parameters['quench_min_temp']))
    if 'quench_step_size' not in parameters:
        parameters['quench_step_size'] = 0.01
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Setting quench_step_size = {0}'.format(parameters['quench_step_size']))
    if 'quench_n_steps_1' not in parameters:
        parameters['quench_n_steps_1'] = 10000
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Setting quench_n_steps_1 = {0}'.format(parameters['quench_n_steps_1']))
    if 'quench_n_steps_2' not in parameters:
        parameters['quench_n_steps_2'] = parameters['quench_n_steps_1']*2
        if rank == 0:
            if 'quench' in parameters['mutation_options']:
                logger.info('Setting quench_n_steps_2 = {0}'.format(parameters['quench_n_steps_2']))
    if 'isolate_mutation' not in parameters:
        parameters['isolate_mutation'] = False
        if rank == 0:
            logger.info('Setting isolate_mutation flag = {0}'.format(parameters['isolate_mutation']))
    
    #Parameters for Selection
    if 'energy_cutoff_factor' not in parameters:
        parameters['energy_cutoff_factor'] = 10.0
        if rank ==0:
            logger.info('Setting energy_cutoff_factor = {0}'.format(parameters['energy_cutoff_factor']))
   
    """
    stem_parameters_keys = ['psf_file','stem_ref','electron_energy','spherical_aberration','defocus','aperture_semiangle','source_size','slice_size','chromatic_aberration_coefficient','delta_E','aber','scale_factor']
    if 'stem' in parameters['fitness_scheme']:
        stemparams = {}
        for one in stem_parameters_keys:
            if one not in parameters :
                logger.critical('Must provide stem_parameters %s for STEM_Cost calculation'%one)
                raise RuntimeError("STEM parameters not specified.  Cannot simulate image files")
            else:
                stemparams[one] = parameters[one] 
        if 'fitting_coeff' not in parameters:
            parameters['fitting_coeff'] = [1]
        else:
            stemparams['fitting_coeff'] = parameters['fitting_coeff']
        if 'grid_sim2exp' not in parameters:
            stemparams['grid_sim2exp'] = 1
        else:
            stemparams['grid_sim2exp'] = parameters['grid_sim2exp']
        if 'pixelshift' not in parameters:
            stemparams['pixelshift'] = False
        else:
            stemparams['pixelshift'] = parameters['pixelshift']
    if 'stem_keep_files' not in parameters:
        parameters['stem_keep_files'] = True
        if rank == 0:
            if 'stem' in parameters['fitness_scheme']:
                logger.info('Setting stem_keep_files = {0}'.format(parameters['stem_keep_files']))
    else:
        parameters['stem_keep_files'] = parameters['stem_keep_files']
    if 'psf_file' in parameters and not parameters['psf_file']==None:
        fileobj = open(parameters['psf_file'],'r')
        lines = fileobj.readlines()
        nk = len(lines)
        parameters['pixels'] = nk
        stemparams['pixels'] = nk
        parameters['psf'] = numpy.empty([nk,nk],dtype=float)
        for x in range(0,nk):
            parameters['psf'][x] = lines[x].split()
        stemparams['psf'] = parameters['psf']
        stemparams['psf_file'] = parameters['psf_file']
    if 'stem_coeff' not in parameters:
        parameters['stem_coeff'] = None
        if rank == 0:
            if 'stem' in parameters['fitness_scheme']:
                logger.info('Setting stem_coeff with first individual')
    else:
        try:
            parameters['stem_coeff'] = float(parameters['stem_coeff'])
        except:
            if rank == 0:
                if 'stem' in parameters['fitness_scheme']:
                    logger.warning('Trouble reading stem_coeff input. stem_coeff = {0}'.format(parameters['stem_coeff']))
            parameters['stem_coeff'] = None
    if 'stem' in parameters['fitness_scheme']:
        #Initialize function for experimental image
        from StructOpt.tools.StemCalc import ConvStem
        logger.info('Initializing ConvStem Calculator')
        parameters['stemcalc'] = ConvStem(parameters=stemparams, tmp_dir='/'+os.getcwd()+'/ConvStemImages/', keep_files=parameters['stem_keep_files'])
    else:
        parameters['stemcalc'] = None
    """
    
    if 'constrain_swaps' not in parameters:
        if 'IntSwap' in parameters['mutation_options']:
            parameters['swaplist'] = None
            if rank == 0:
                logger.info('Setting swaplist for IntSwap = None')
        else:
            parameters['swaplist'] = False
            parameters['constrain_swaps'] = False
            if rank == 0:
                logger.info('Setting swaplist = False')
    if 'natural_selection_scheme' not in parameters:
        parameters['natural_selection_scheme'] = natselectscheme
        if rank ==0:
            logger.info('Setting natural_selection_scheme = {0}'.format(parameters['natural_selection_scheme']))
    if 'tournsize' not in parameters:
        parameters['tournsize'] = 3
        if 'tournament' in parameters['selection_scheme'] or 'tournament' in parameters['natural_selection_scheme']:
            if rank ==0:
                logger.info('Setting Tournament size (tournsize) = {0}'.format(parameters['tournsize']))
    if 'fusslimit' not in parameters:
        parameters['fusslimit'] = 10.0
        if rank ==0:
            logger.info('Setting FUSS limit (fusslimit) = {0}'.format(parameters['fusslimit']))
    if 'metropolis_temp' not in parameters:
        parameters['metropolis_temp'] = 30.0
        if rank == 0:
            logger.info('Setting metropolis_temp = {0}'.format(parameters['metropolis_temp']))
    if 'mark' not in parameters:
        parameters['mark'] = None
    #Parameters for Convergence
    if 'convergence_scheme' not in parameters:
        parameters['convergence_scheme'] = 'max_gen'
        if rank ==0:
            logger.info('Setting convergence scheme (convergence_scheme) = {0}'.format(parameters['convergence_scheme']))
    if 'maxgen' not in parameters:
        parameters['maxgen'] = 5
        if rank ==0:
            logger.info('Setting Max Number of generations (maxgen) = {0}'.format(parameters['maxgen']))
    if 'reqrep' not in parameters:
        parameters['reqrep'] = 10
        if rank ==0:
            if 'rep' in parameters['convergence_scheme']:
                logger.info('Setting max number of energy repetitions (reqrep) = {0}'.format(parameters['reqrep']))
    if 'tolerance' not in parameters:
        parameters['tolerance'] = 0.001
        if rank ==0:
            if 'rep' in parameters['convergence_scheme']:
                logger.info('Setting energy tolerance (tolerance) = {0}'.format(parameters['tolerance']))
    if 'predator' not in parameters:
        parameters['predator'] = 'mutation_dups'
        if rank == 0:
            logger.info('Setting predator = {0}'.format(parameters['predator']))
    if 'adaptbegin' not in parameters:
        parameters['adaptbegin'] = 0.75
        if rank == 0:
            if parameters['predator'] == 'adapting':
                logger.info('Setting adaptation predator to begin (adaptbegin) at genrep*{0}'.format(parameters['adaptbegin']))
    if 'adaptmultiplier' not in parameters:
        parameters['adaptmultiplier'] = 3.0
        if rank == 0:
            if parameters['predator'] == 'adapting':
                logger.info('Setting adaptation predator multiplier (adaptmultiplier) = {0}'.format(parameters['adaptmulitplier']))
    if 'demin' not in parameters:
        parameters['demin'] = 0.005
        if rank == 0:
            logger.info('Setting cutoff convergence energy (demin) = {0}'.format(parameters['demin']))
    
    return parameters

class dummy_logger():
    def __init__(self):
        return
    def critical(self,message):
        print 'CRITICAL: {0}'.format(message)
        return
    def debug(self,message):
        print 'DEBUG: {0}'.format(message)
        return
    def warning(self,message):
        print 'WARNING: {0}'.format(message)
        return
    def info(self,message):
        print 'MESSAGE: {0}'.format(message)
        return

class dummy_logger_no_write():
    def __init__(self):
        return
    def critical(self,message):
        return
    def debug(self,message):
        return
    def warning(self,message):
        return
    def info(self,message):
        return
