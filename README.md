### Description

This [snakemake](https://snakemake.readthedocs.io/en/stable/) workflow
implements a fitting pipeline for the NTDMC trachoma model. Starting
from geostatistical infection data, this workflow simulates the spread
of trachoma for an (statistical) ensemble of models with different
beta parameter values, for each implementation unit present in the
original data.

The simplified workflow can visualised as:

![simplified dag](./img/simplified_dag.svg)

The workflow can be divided in two successive parts.

1. The AMIS algorithm is used to generate an ensemble of values for
the transmission parameter (beta), along with their corresponding
statistical weight. In practice, the workflow splits the IUs into
groups processed in parallel.  IUs are grouped together according to
values of first MDA year, last MDA year, mean prevalence level. For
instance, IU subset (`first_mda`=2008, `last_mda`=2019, `level`=3)
groups together IUs for which the first MDA occurs in year 2008, the
last MDA in year 2019, and for which the mean infection prevalence is
between 20% and 30% (assuming a prevalence level spans a 10%
prevalence range).

2. The epidemiological model is simulated for a subset of the
generated parameter values, sampled according to their weight. For
each implementation unit, the model is first simulated over past years
up to an arbitraty end year -- typically the present year. The model
is then forward simulated from this end year to an arbitrary point in
time in the future.

A more accurate visualisation of the workflow is therefore, with 3
implementation units belonging to two different groups

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

To execute the pipeline for all implementation units in the input data, run
```shell
snakemake --cores N 
```

replacing `N` by the number of cores available. As much as possible,
depending on the number of available cores, subsets of IUs will be
processed in parallel.

It is possible to process individual implementation units, for
instance

```shell
snakemake --cores 4 results/infection_ETH18604.cdv
```

Or several of them (using bash shell brace expansion notation)
```shell
snakemake --cores 4 data/infection_{ETH18604,ETH18551,ETH18644}.csv
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
`{FIRST_MDA}`, `{LAST_MDA}`, `{GROUP}` or `{IUCODE}`. For example:

```python
def get_input_amis_file(wildcards):
    start_MDA, last_MDA, group = get_group_from_IU_code(wildcards)
    return f"results/amis_output_{start_MDA}_{last_MDA}_{group}.csv"


rule sample_parameter_values:
    input:
        get_input_amis_file,
    output:
        "results/sampled_parameters_{IUCODE}.csv",
    params:
        nsamples=config["nsamples_beta"],
    script:
        "scripts/sample_parameters.py"
```

The above rule tells snakemake how to produce files named as
`sampled_parameters_{IUCODE}.csv`, `{IUCODE}` being a placeholder for
an arbitraty string. For example,

```shell
snakemake results/sampled_parameters_ETH18644.csv
```

will execute rule `sammple_parameter_values` with `IUCODE=ETH18644`.

Moreover, input files are defined as the return value of input
function `get_input_amis_file`, see [Input
functions](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#input-functions). This
function is passed a `wildcard` object defined by Snakemake, carrying
the value of wildcard `IUCODE` as an attribute
(`wildcards.IUCODE`). From this `wildcard` object, function
`get_group_from_IU_code` is able to return the corresponding IU group
this particular IU belongs to. Note that the evaluation of input
functions is done when evaluating the graph of jobs.

#### Target rule `all`

Running 

```shell
snakemake
```

will execute the target (default) rule `all`. It only consists of an
input function aggregating model outputs for all implementation units
in the data.

#### Checkpoint rule

Rule `group_ius` *must* have been executed for function
`get_subset_from_IU_code` to be able to assign individual
implemenation unit codes to their respecitve group. This rule is
therefore defined as a checkpoint rule, see [Data-dependent
conditional
execution](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?highlight=checkpoint#data-dependent-conditional-execution). In
practice, Snakemake will register that rules depending on the
execution of function `get_subset_from_IU_code` require the output of
rule `group_ius`, execute `group_ius`, and re-evaluate the dependency
graph.
