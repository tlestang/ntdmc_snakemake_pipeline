.. _required-software:

Required software
=================

Snakemake (required)
    To execute the workflow.
Git (recommended)
    To download and update workflow definition and support files.
Conda (recommended)
    To enable automatic creation of computational environments.

The Snakemake workflow directly depends on

- Pandas
- The NTDMC trachoma python model
- The AMIS R package.

Using snakemake together with the conda package (see ??) manager
allows the automatic installation of the dependencies listed above,
and is therefore the recommended approach for most users.
