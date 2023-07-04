Quickstart
==========

Starting dataset
----------------

Our starting point is statistical data about the prevalence of
trachoma at several locations across Ethiopia, Nigeria, Uganda,
Burkina Faso, Tanzania, Cameroon and Senegal.  Each location is
refered to as an *Implementation Unit (IU)* and is uniquely described
by a IU code (e.g. ETH18538 or SEN40173).  Prevalence statistics are
given after several rounds of *Mass Drug Application*, described in
the starting data file.  In this case, MDAs are carried yearly between
values indicated in columns ``start_MDA`` and ``last_MDA``.

Configuration
--------------

Pipeline parameters can be set by editing the ``config.yaml`` file.
For the purpose of this fast-paced introduction the only parameters we
are concerned with are

data
    Location of the dataset, relative to the workflow file (``Snakefile``).
end_sim_year
    Model simulation end year.
nsamples_beta
    Number of parameter values to simulate.

Executing the workflow
----------------------

The pipeline generates $N$ output from NTDMC trachoma model (disease
infection and prevalence) for an ensemble of $N$ epidemiological
parameter values, for each IU in the dataset, simulated up to an
arbitrary (user-defined) final year.  The number of simulated
parameter values per IU is also user defined.

Simulation output can be generated for a single IU or a set of IUs:

Single IU:

.. code-block:: shell

   snakemake --cores 1 --use-conda --conda-frontend conda \
     results/infection_SEN40173.csv

Multiple IUs

.. code-block:: shell

   snakemake --cores 1 --use-conda --conda-frontend \
     results/infection_{ETH18587,ETH18542,TZA46580}.csv

Alternatively, simulationp output can generated for all IUs in the
dataset with:

.. code-block:: shell

   snakemake --cores 1 --use-conda --conda-frontend conda

.. note::

   The ``--use-conda`` and ``--conda-frontend`` options enable the
   automatic creation and activation of conda environments. See
   `(Snakemake)Integrated Package Management
   <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management>`_.
   Should you choose to manage the installation the required workflow
   dependencies yourself (see :ref:`required-software`), both options
   can be ommitted.

Partial execution
-----------------

The workflow generates a number of intermediate results.  For
instance, you can ask snakemake to generate the sampled parameter
values.

snakemake --cores 1 --use-conda --conda-frontend \
     results/sampled_parameters_{ETH18587,ETH18542,TZA46580}.csv

In this case, the workflow execution will stop at the generation of
the parameter values, and will not carry on with the model
simulation. See :ref:`Worflow results` for a list of all available
workflow results.


Specifying workflow parameters at the command line
--------------------------------------------------

Parameter values specified in the configuration file (``config.yaml``)
can be overriden using the ``--config`` option:

.. code:: shell

   snakemake --cores 1 --use-conda --conda-frontend conda



