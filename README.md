### Description

This [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow
implements a fitting pipeline for the NTDMC trachoma model. Starting
from geostatistical infection data, this workflow simul spread of
trachoma for an (statistical) ensemble of models with different beta
parameter values, for each implementation unit present in the original
data.

The simplified workflow can visualised as:

![simplified dag](./img/simplified_dag.svg)

In practice, the workflow splits the IUs into groups processed in
parallel.  IUs are grouped together according to values of first MDA
year, last MDA year, mean prevalence level. For instance, IU subset
(`first_mda`=2008, `last_mda`=2019, `level`=3) groups together IUs for
which the first MDA occurs in year 2008, the last MDA in year 2019,
and for which the mean infection prevalence is between 20% and 30%
(assuming a prevalence level spans a 10% prevalence range).

A more accurate visualisation of the workflow is therefore, with 3 IU
groups:

![parallel_dag](./img/parallel_dag.svg)

### Installation

You'll need 
- A Python 3.8+ installation with the `scipy` and `pandas` packages availabe.
- A R installation with the [trachomAMIS
  package](https://github.com/OxfordRSE/trachomAMIS) available.

For instance, you could install the miniconda python distribution and
create a `trachoma` environment:
```shell
(base) $ conda create --name trachoma scipy pandas
```
Then, install snakemake in the `trachoma` environment
```shell
(base) $ conda activate trachoma
(trachoma) $ pip install snakemake
```

Finally, install the trachomAMIS R package

```R
devtools::install_github("OxfordRSE/trachomAMIS")
```

### Usage

To execute the pipeline for all IU subsets in the input data, run
```shell
snakemake --cores N 
```

replacing `N` by the number of cores available. As much as possible,
depending on the number of available cores, subsets of IUs will be
processed in parallel.

It is possible to process individual subsets, for instance

```shell
snakemake --cores 4 data/model_output_2008_2019_level_3
```

Or several of them (using bash shell brace expansion notation)
```shell
snakemake --cores 4 data/model_output_{2008_2019_level_3,2008_2017_level_2}
```
### Under the hood

The fitting pipeline is defined in the [Snakefile](./Snakefile). This
file breaks down the pipeline into a set of rules with input data and
output data. Snakemake processes the Snakefile and builds a dependency
graph connecting rules that procude output that other rules consume as
input. Snakemake then works this graph backward, executing rules in
the right order to produce the final output.

#### Wildcards

Most of the rules in the workflow are based on *wildcards*, such as
`{FIRST_MDA}`, `{LAST_MDA}` or `{GROUP}`. For example:

```python
rule sample_parameter_values:
    input:
        "amis_output_{FIRST_MDA}_{LAST_MDA}_{GROUP}.csv",
    output:
        "sampled_parameters_{FIRST_MDA}_{LAST_MDA}_{GROUP}.csv",
    params:
        nsamples=10,
    script:
        "scripts/sample_parameters.py"
```

The above rule tells snakemake how to produce files named as
`sampled_parameters_{FIRST_MDA}_{LAST_MDA}_{GROUP}.csv`, given input
file `amis_output_{FIRST_MDA}_{LAST_MDA}_{GROUP}.csv`, for any value
of the three wildcards. For instance,

```shell
snakemake sampled_parameters_2008_2018_2.csv
```

will try to execute rule `sample_parameter_values` to produce
`sampled_parameters_2008_2018_2.csv` from
`amis_output_2008_2018_2.csv`, if it exists.

#### The target rule `all`

Running 

```shell
snakemake
```

will execute the target (default) rule `all`.

```python
def aggregate_input(wildcards):
    from pandas import read_csv

    checkpoint_output = checkpoints.group_ius.get(**wildcards).output[0]
    grouped = read_csv(checkpoint_output).groupby(["start_MDA", "last_MDA", "group"])
    return [f"data/model_output_{name[0]}_{name[1]}_{name[2]}" for name, _ in grouped]


rule all:
    input:
        aggregate_input,
```

Its input is defined as a Python function that builds and return the
list of model output directories for each group of IUs (`first_mda`,
`last_mda`, `level`). It does so by reading the input data augmented
with the `level` column (the output of the `group_ius` rule) and
grouping rows by values of columns `first_mda`, `last_mda` and
`level`.

Because, the output of rule `group_ius` *must* have been created to
determine input of ruel `all`, rule `group_ius` is defined as a
checkpoint rule. In practice, Snakemake will register that rule `all`
requires the output of rule `group_ius`, execute `group_ius`, and
re-evaluate the dependency graph. See [Data-dependent conditional
execution](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?highlight=checkpoint#data-dependent-conditional-execution)
in the Snakemake docs.
