==========
JobManager
==========




.. _sec-introduction:

Introduction
------------

The purpose of the ``JobManager`` module is to provide a python wrapper for submitting and tracking jobs in a queue environment.

.. _sec-configuration:

Configuration
-------------

The ``JobManager`` is initially built for a PBS queue environment, so many of the commands will have to be modified for usage in a different queue environment. These customizations will likely take place in the following files.

1. The ``submit`` and ``write_submit`` function in the ``structopt/utilities/job_manager.py`` file will likely need to be updated to reflect your specific queue environment.

2. The dictionaries held in ``structopt.utilities/rc.py`` is the first attempt to store some dictionaries specific to the queue environment. Many queue specific variables are drawn from here.

.. _sec-submit:

Submitting jobs
---------------

.. _sec-submit-single:

Single job
~~~~~~~~~~

The script below is an example script of submitting a single job to a queue using the ``JobManager``. The optimization run is a short run of a Au\ :sub:`55`\ nanoparticle using only LAMMPS. A large part of the script is defining the input, which goes into the ``JobManager`` class. These inputs are given below.

1. ``calcdir``: This is a string that tells where the calculation is run. Note that the calculation itself is run within the ``calcdir/logs{time}`` directory, which is created when the job starts to run on the queue. Unless an absolute path, the ``calcdir`` directory is always given with respect to directory that the job script is run from

2. ``optimizer``: This is a string of the optimizer file used for the calculation. These files can be found in the ``structopt/optimizers`` folder. Upon run, a copy of this script is placed insde of the ``calcdir`` directory and accessed from there.

3. ``structopt_parameters``: This is a dictionary object that should mirror the input file you are trying to submit

4. ``submit_parameters``: This dictionary holds the submit parameters. These will be specific to the queue system in use. In this example, we specify the the submission system, queue, number of nodes, number of cores, and walltime.

.. code-block:: python

    from structopt.utilities.job_manager import JobManager
    from structopt.utilities.exceptions import Running, Submitted, Queued

    calcdir = 'job_manager_examples/Au55-example'

    structopt_parameters = {
        "seed": 0,
        "structure_type": "cluster",
        ...
    }

    submit_parameters = {'system': 'PBS',
                         'queue': 'morgan2',
                         'nodes': 1,
                         'cores': 12,
                         'walltime': 12}

    optimizer = 'genetic.py'

    job = JobManager(calcdir, optimizer, structopt_parameters, submit_parameters)
    job.optimize()

Upon running this script, the user should get back an exception called ``structopt.utilities.exceptions.Submitted`` with the jobid. This is normal behavior and communicates that the job has successfully been submitted.

.. _sec-submit-multiple:

Multiple jobs
~~~~~~~~~~~~~

One advantage of the job manager is that it allows one to submit multiple jobs to the queue. This is often useful for tuning the optimizer against different inputs. The script below is an example of submitting the same job at different seeds.

In the previous script, submitting a single job successfully with ``JobManager.optimizer`` method resulted in an exception. We can catch this exception with a ``try`` and ``except`` statement. This is shown below in the script where upon a successful submission, the script prints out the jobid to the user.

.. code-block:: python

    from structopt.utilities.job_manager import JobManager
    from structopt.utilities.exceptions import Running, Submitted, Queued

    structopt_parameters = {
        "seed": 0,
        "structure_type": "cluster",
        ...
    }

    submit_parameters = {'system': 'PBS',
                         'queue': 'morgan2',
                         'nodes': 1,
                         'cores': 12,
                         'walltime': 12}

    optimizer = 'genetic.py'

    seeds = [0, 1, 2, 3, 4]
    for seed in seeds:
        structopt_parameters['seed'] = seed
        calcdir = 'job_manager_examples/Au55-seed-{}'.format(seed)

        job = JobManager(calcdir, optimizer, structopt_parameters, submit_parameters)

        try:
            job.optimize()
        except Submitted:
            print(calcdir, job.get_jobid(), 'submitted')

::

    job_manager_examples/Au55-seed-0 936454.bardeen.msae.wisc.edu submitted
    job_manager_examples/Au55-seed-1 936455.bardeen.msae.wisc.edu submitted
    job_manager_examples/Au55-seed-2 936456.bardeen.msae.wisc.edu submitted
    job_manager_examples/Au55-seed-3 936457.bardeen.msae.wisc.edu submitted
    job_manager_examples/Au55-seed-4 936458.bardeen.msae.wisc.edu submitted

.. _sec-track:

Tracking jobs
-------------

In the previous section, we covered how to submit a new job from an empty directory. This is done by first initializing an instance of the ``StructOpt.utilities.job_manager.JobManager`` class with a calculation directory along with some input files and then submitting the job with the ``JobManager.optimize`` method. The ``JobManager.optimize`` method knows what to do because upon initialization, it detected an empty directory. If the directory was not empty and contained a StructOpt job, the ``JobManager`` knows what to do with it if ``optimize`` was run again. This is all done with exceptions.

The four primary exceptions that are returned upon executing the ``optimize`` method are below along with their explanations.

1. ``Submitted``: This exception is returned if a job is submitted from the directory. This is done when ``JobManager.optimize`` is called in an empty directory or ``JobManager.optimize`` is called with the kwarg ``restart=True`` in a directory where a job is not queued or running.

2. ``Queued``: The job is queued and has not started running. There should be no output files to be analyzed.

3. ``Running``: The job is running and output files should be continously be updated. These output files can be used for analysis before the job has finished running.

4. ``UnknownState``: This is returned if the ``calcdir`` is not an empty directory doesn't detect it as a StructOpt run. A StructOpt run is detected when a ``structopt.in.json`` file is found in the ``calcdir``.

Note that if no exception is returned, it means the job is done and is ready to be analyzed. ``Job.optimize`` does nothing in this case.

One way of using these three exceptions is below. If the job is submitted or Queued, we want the script to stop and not submit the job. If it is running, additional commands can be used to track the progress of the job.

.. code-block:: python

    from structopt.utilities.job_manager import JobManager
    from structopt.utilities.exceptions import Running, Submitted, Queued

    calcdir = 'job_manager_examples/Au55-example'

    structopt_parameters = {
        "seed": 0,
        "structure_type": "cluster",
        ...
    }

    submit_parameters = {'system': 'PBS',
                         'queue': 'morgan2',
                         'nodes': 1,
                         'cores': 12,
                         'walltime': 12}

    optimizer = 'genetic.py'

    job = JobManager(calcdir, optimizer, structopt_parameters, submit_parameters)
    try:
        job.optimize()
    except (Submitted, Queued):
        print(calcdir, job.get_jobid(), 'submitted or queued')
    except Running:
        pass

::

    job_manager_examples/Au55-example 936453.bardeen.msae.wisc.edu submitted or queued

.. _sec-restart:

Restarting jobs
---------------

Sometimes jobs need to be restarted or continued from the last generation. The ``JobManager`` does this by submitting a new job from the same ``calcdir`` folder the previous job was run in. Because calculations take place in unique ``log{time}`` directories, the job will run in a new ``log{time}`` directory. Furthermore, the ``JobManager`` modifies the ``structopt.in.json`` file so the initial population of the new job are the XYZ files of the last generation of the previous run.  The code below is an example of restarting the first run of this example. The only difference between this code and the one presented in the previous section is that a ``restart=True`` kwarg has been added to the ``JobManager.optimize`` command.

.. code-block:: python

    from structopt.utilities.job_manager import JobManager
    from structopt.utilities.exceptions import Running, Submitted, Queued

    calcdir = 'job_manager_examples/Au55-example'

    structopt_parameters = {
        "seed": 0,
        "structure_type": "aperiodic",
        ...
    }

    submit_parameters = {'system': 'PBS',
                         'queue': 'morgan2',
                         'nodes': 1,
                         'cores': 12,
                         'walltime': 12}

    optimizer = 'genetic.py'

    job = JobManager(calcdir, optimizer, structopt_parameters, submit_parameters)
    job.optimize(restart=True)
