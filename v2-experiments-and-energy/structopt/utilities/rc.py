QUEUE_OPTIONS = {'PBS': {'prefix': '#PBS',
                         'submit': 'qsub',
                         'job_name': '-N {}',
                         'queue': '-q {}',
                         'nodes_cores': '-l nodes={0}:ppn={1}',
                         'walltime': '-l walltime={}:00:00',
                         'misc': '-j oe'}}

RUN_OPTIONS = {'mpirun': 'mpirun',
               'python': 'python'}

CUSTOM_LINES = '''source activate py35
export PYTHONPATH=$HOME/research/StructOpt_modular/:$PYTHONPATH
export PATH=/share/apps/openmpi-1.10.0_no_ib/bin/:$PATH\n'''

