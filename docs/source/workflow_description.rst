Workflow description
====================

The worflow can be divided into two distinct stages:

- Parameter inference
- Model simulation

The parameter inference stage is about estimating a probablity
distribution for values of the epidemiological model parameters, using
the Adaptive Multilevel Importance Sampling algorithm. This stage
yields a collection of parameter values, together with their
corresponding statistical weight.

The model simulation stage then forward simulates the model for a
arbitrary number of parameter values sampled in the ensemble generated
by the inference stage.

the AMIS
algorithm

Parameter inference stage
-------------------------

IU grouping
###########

The parameter inference (AMIS) algorithm is executed over several
groups of IUs, instead of the whole dataset.  This allows to divide
the computation into independent pieces of work that can run in
parallel.

Currently, IUs are grouped together according to values of first MDA
year, last MDA year, mean prevalence level. For instance, IU a group
``(first_mda=2008, last_mda=2019, level=3)1`` groups together IUs for
which the first MDA occurs in year 2008, the last MDA in year 2019,
and for which the mean infection prevalence is between 20% and 30%
(assuming a prevalence level spans a 10% prevalence range).  IU
grouping will likely become configurable in the near future.

Prevalence maps
###############

The input of the AMIS algorithm is a **prevalence map**, which is
derived from the original dataset.  A prevalence map contains K
samples of prevalence value for each IU in the dataset (or IU
group). It is stored in the CSV format with N rows and K columns.
