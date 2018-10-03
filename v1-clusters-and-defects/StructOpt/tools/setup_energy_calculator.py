from StructOpt.tools.lammps import LAMMPS
import os
import json
try:
    from mpi4py import MPI
except ImportError:
    pass
import logging

def setup_energy_calculator(Optimizer,mod,relax):
    if mod == 'VASP':
        args = json.load(open('vasp_inp.json'))
        return
    elif mod == 'LAMMPS':
        if 'SetCalc' in Optimizer.debug:
            debug=True
            logger = logging.getLogger(Optimizer.loggername)
        else:
            debug=False
        
        args = json.load(open('lammps_inp.json'))
        atomlist=Optimizer.atomlist
        atomlist=sorted(atomlist,key=lambda symbol: symbol[0])
        if args["pair_style"]=='tersoff':
            if debug:
                logger.info('Setting up LAMMPS calculator with Tersoff potential')
            parcoff = '* * {0}'.format(args["pot_file"])
            for one in atomlist:
                parcoff+=' {0}'.format(one[0])
            pair_coeff = [parcoff]
            mass = ['1 {0}'.format(atomlist[0][2])]
            if len(atomlist) > 1:
                for i in range(len(atomlist)-1):
                    mass.append('{0} {1}'.format(i+2,atomlist[i+1][2]))
            parameters = { 'pair_style' : args["pair_style"], \
            'pair_coeff' : pair_coeff , 'mass' : mass }
            filesL = [ args["pot_file"] ]
        elif args["pair_style"] =='eam':
            if debug:
                logger.info('Setting up LAMMPS calculator with EAM potential')
            pair_coeff = [ '* * {0}'.format(args["pot_file"])]
            parameters = { 'pair_style' : args["pair_style"], 'pair_coeff' : pair_coeff }
            filesL = [ args["pot_file"] ]
        elif args["pair_style"]=='eam/fs':
            if debug:
                logger.info('Setting up LAMMPS calculator with EAM/FS potential')
            parcoff = '* * {0}'.format(args["pot_file"])
            for one in atomlist:
                parcoff+=' {0}'.format(one[0])
            pair_coeff = [parcoff]
            mass = ['1 {0}'.format(atomlist[0][2])]
            if len(atomlist) > 1:
                for i in range(len(atomlist)-1):
                    mass.append('{0} {1}'.format(i+2,atomlist[i+1][2]))
            parameters = { 'pair_style' : args["pair_style"], 
            'pair_coeff' : pair_coeff , 'mass' : mass }
            filesL = [ args["pot_file"] ]
        elif args["pair_style"]=='eam/alloy':
            if debug:
                logger.info('Setting up LAMMPS calculator with EAM/ALLOY potential')
            parcoff = '* * {0}'.format(args["pot_file"])
            for one in atomlist:
                parcoff+=' {0}'.format(one[0])
            pair_coeff = [parcoff]
            parameters = { 'pair_style' : args["pair_style"], 'pair_coeff' : pair_coeff}
            filesL = [ args["pot_file"] ]
        elif args["pair_style"]=='eam/cd':
            if debug:
                logger.info('Setting up LAMMPS calculator with EAM/CD potential')
            parcoff = '* * {0}'.format(args["pot_file"])
            for one in atomlist:
                parcoff+=' {0}'.format(one[0])
            pair_coeff = [parcoff]
            parameters = { 'pair_style' : args["pair_style"], 'pair_coeff' : pair_coeff}
            filesL = [ args["pot_file"] ]
        elif args["pair_style"]=='edip':
            if debug:
                logger.info('Setting up LAMMPS calculator with EDIP potential')
            parcoff = '* * {0}'.format(args["pot_file"])
            for one in atomlist:
                parcoff+=' {0}'.format(one[0])
            pair_coeff = [parcoff]
            mass = ['1 {0}'.format(atomlist[0][2])]
            if len(atomlist) > 1:
                for i in range(len(atomlist)-1):
                    mass.append('{0} {1}'.format(i+2,atomlist[i+1][2]))
            parameters = { 'pair_style' : args["pair_style"], \
            'pair_coeff' : pair_coeff , 'mass' : mass, 'newton': 'on' }
            filesL = [ args["pot_file"] ]
        elif args["pair_style"]=='bop':
            if debug:
                logger.info('Setting up LAMMPS calculator with BOP potential')
            parcoff = '* * {0}'.format(args["pot_file"])
            for one in atomlist:
                parcoff+=' {0}'.format(one[0])
            parcoff+='\ncommunicate single cutoff {0}'.format(Optimizer.bopcutoff)
            pair_coeff = [parcoff]
            mass = ['1 {0}'.format(atomlist[0][2])]
            if len(atomlist) > 1:
                for i in range(len(atomlist)-1):
                    mass.append('{0} {1}'.format(i+2,atomlist[i+1][2]))
            parameters = { 'pair_style' : args["pair_style"], \
            'pair_coeff' : pair_coeff, 'mass' : mass, 'newton': 'on' }
            filesL = [ args["pot_file"] ]
        elif args["pair_style"]=='buck':
            if debug:
                logger.info('Setting up LAMMPS calculator with Buckingham potential')
            pairstyle='{0} {1}'.format(args["pair_style"], Optimizer.buckcutoff)
            pair_coeff=Optimizer.buckparameters
            mass = ['1 {0}'.format(atomlist[0][2])]
            if len(atomlist) > 1:
                for i in range(len(atomlist)-1):
                    mass.append('{0} {1}'.format(i+2, atomlist[i+1][2]))    
            parameters = {'pair_style': pairstyle, 'pair_coeff': pair_coeff, \
             'mass' : mass }
            filesL=None
        ## ZS
        elif 'lj' in args["pair_style"]:
            if debug:
                logger.info('Setting up LAMMPS calculator with Lennard Jones potential')
            pairstyle=args["pair_style"]
            pair_coeff=[args["pair_coeff"]]
            try:
                cutoff = float(pairstyle.split()[-1])
            except ValueError:
                if len(pair_coeff[0].split())==4:
                    cutoff = pair_coeff[0].split()[3]*5
                elif len(pair_coeff[0].split())==5:
                    cutoff = pair_coeff[0].split()[4]
                pairstyle = "%s %s"%(pairstyle,cutoff)
            mass = ['1 {0}'.format(atomlist[0][2])]
            if len(atomlist) > 1:
                for i in range(len(atomlist)-1):
                    mass.append('{0} {1}'.format(i+2,atomlist[i+1][2]))
            parameters = { 'pair_style' : pairstyle, \
            'pair_coeff' : pair_coeff, 'mass' : mass}
            filesL=None
        
        elif args["pair_style"]=='other':
            """WARNING: This style still needs work. Intended to allow user flexibility with potential specification"""
            if debug:
                logger.info('Setting up LAMMPS calculator with user input potential')
            mass = ['1 {0}'.format(atomlist[0][2])]
            if len(atomlist) > 1:
                for i in range(len(atomlist)-1):
                    mass.append('{0} {1}'.format(i+2,atomlist[i+1][2]))
            if args["ps_other"]!=None:
                if 'newton' in args["ps_other"]:
                    parameters = {'pair_style' : args["ps_name"], \
                    'pair_coeff': [args["pair_coeff"]], 'mass': mass,'newton':'on'}
                else:
                    parameters = {'pair_style' : args["ps_name"], \
                    'pair_coeff': [args["pair_coeff"]], 'mass': mass,'newton':'off'}
                if 'charges' in args["ps_other"]:
                    cs=args["ps_other"].split('charges:')
                    parameters['mass'][len(parameters['mass'])-1] +=cs[1]
                    parameters['newton']+='\natom_style charge'
            else:
                parameters = {'pair_style' : args["ps_name"], \
                'pair_coeff': args["pair_coeff"], 'mass': mass}
            if args["pot_file"] !=None:
                filesL = [ args["pot_file"] ]
            else:
                filesL=None
        else:
            if debug:
                logger.warn('No LAMMPS potential recognized. Setting up LAMMPS calculator with Lennard Jones potential')
            parameters={}
            filesL=None
            print 'WARNING: No LAMMPS potential recognized. Assuming Lennard Jones Potential'
        if args["minimize"] != None:
            if debug:
                logger.info('Adding local energy minimizer to LAMMPS calculator')
            try:
                parameters['mass'][len(parameters['mass'])-1] += '\nmin_style {0}'.format(args["min_style"])
            except KeyError:
                parameters['pair_coeff'][0] += '\nmin_style {0}'.format(args["min_style"])
            parameters['minimize'] = args["minimize"]
        
        if not relax:
            parameters['minimize'] = "1e-8 1e-8 0 0"
        parameters['thermosteps'] = args["thermo_steps"]
        if args["keep_files"]:
            if debug:
                logger.info('Setting up directory for keeping LAMMPS files')
            try:
                rank = MPI.COMM_WORLD.Get_rank()
                if Optimizer.parallel:
                    if 'Island_Method' not in Optimizer.algorithm_type:
                        real_rank = MPI.COMM_WORLD.Get_rank()
                        rank = 0
                        path = os.path.join(os.getcwd(),'{0}-rank{1}'.format(Optimizer.filename,rank))
                        if not os.path.exists(os.path.join(path,'LAMMPSFiles')):
                            os.mkdir(os.path.join(path,'LAMMPSFiles'))
                            logger.info('Making directory: {0}'.format(os.path.join(path,'LAMMPSFiles')))
            except:
                rank = 0
            if filesL != None:
                path = os.path.join(os.getcwd(),'{0}-rank{1}'.format(Optimizer.filename,rank))
                if Optimizer.parallel and ('Island_Method' not in Optimizer.algorithm_type):
                    tmpdir = os.path.join(os.path.join(path, 'LAMMPSFiles'),'rank-{0}'.format(real_rank))
                    calc = LAMMPS(parameters=parameters, files=filesL, \
                        keep_tmp_files=True, tmp_dir=tmpdir)
                else:
                    calc = LAMMPS(parameters=parameters, files=filesL, \
                        keep_tmp_files=True, tmp_dir=os.path.join(path, 'LAMMPSFiles'))
            else:
                path = os.path.join(os.getcwd(),'{0}-rank{1}'.format(Optimizer.filename,rank))
                #calc = LAMMPS(parameters=parameters, keep_tmp_files=True, \
                #    tmp_dir=os.path.join(path,'LAMMPSFiles'))
                if Optimizer.parallel and ('Island_Method' not in Optimizer.algorithm_type):
                    tmpdir = os.path.join(os.path.join(path, 'LAMMPSFiles'),'rank-{0}'.format(real_rank))
                    calc = LAMMPS(parameters=parameters,  \
                        keep_tmp_files=True, tmp_dir=tmpdir)
                else:
                    calc = LAMMPS(parameters=parameters,  \
                        keep_tmp_files=True, tmp_dir=os.path.join(path, 'LAMMPSFiles'))
        else:
            if filesL != None: 
                calc = LAMMPS(parameters=parameters, files=filesL)
            else:
                calc = LAMMPS(parameters=parameters)
        return calc
