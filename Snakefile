configfile: "config.yaml"


def aggregate_input(wildcards):
    from pandas import read_csv

    checkpoint_output = checkpoints.group_ius.get(**wildcards).output[0]
    grouped = read_csv(checkpoint_output).groupby(["start_MDA", "last_MDA", "group"])
    return [f"data/model_output_{name[0]}_{name[1]}_{name[2]}" for name, _ in grouped]


rule all:
    input:
        aggregate_input,


checkpoint group_ius:
    """
    Group each iu in input data together according to mean prevalence value.
    For example
    group 1: 0 < prev <= 10%
    group 2: 10 < prev <= 20%
    ...
    input:
        CSV data describing one IU per row, and columns 'Logit' and 'Sds'
        describing logit and standard deviation values for each IU.
    output:
        Input CSV data augmented with extra column 'group'.
    """
    input:
        "data/FinalDataTest.csv",
    output:
        "data/FinalDataGroup.csv",
    script:
        "scripts/group_ius.py"


rule make_prevalence_maps:
    """
    Sample params['nsamples'] prevalence values for each IU in
    a (first_mda, last_mda, group) group.
    input:
        CSV data describing one IU per row:

        Logit,Sds,start_MDA,last_MDA,group
        -1.814302776,0.243590266,2008,2017,group_2
        ...

    output:
        CSV data with one row per IU in (first_mda, last_mda, group) group
        and params["nsamples"] columns.
    """
    input:
        "data/FinalDataGroup.csv",
    output:
        "data/prev_map_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv",
    params:
        nsamples=config["nsamples_prevalence_map"],
    script:
        "scripts/compute_prevalence_maps.py"


rule make_mda_file:
    """
    Prepare input file for trachoma model for a specific (first_mda, last_mda)
    subset of IUs. Describes simulation start year, as well as MDA boundary years.
    Simulation end year is described by params["end_sim_year"].
    input:
        CSV data describing one IU per row
    output:
        CSV data

        start_sim_year,end_sim_year,first_mda,last_mda
        2008,2019,2008,2017
    """
    input:
        "data/FinalDataTest.csv",
    output:
        "data/mda_input_{FIRST_MDA}_{LAST_MDA}.csv",
    params:
        end_sim_year=config["end_sim_year"],
    script:
        "scripts/make_mda_files.py"


rule estimate_parameter_weights:
    """
    Runs AMIS for a given (first_mda, last_mda, group) subset of IUs.
    Input:
        prevalence map: CSV data with N prevalence value samples
        per IU, one row per IU.
        mda_input: CSV data describing simulation start and end year,
        as well as MDA first and last year.
    Output:
        CSV data with columns:
            - seeds: Initial seed for trachoma model
            - beta: Sampled beta parameter value
            - sim_prev: Simulated end_year prevalence value
            - then one column per IU in (fisrt_mda, last_mda, group)
              subset, values being statistical weight of beta parameter.
        One row per AMIS iteration.
    """
    input:
        "data/mda_input_{FIRST_MDA}_{LAST_MDA}.csv",
        "data/prev_map_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv",
    output:
        "data/amis_output_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv",
    params:
        nsamples=config["AMIS"]["nsamples"],
        delta=config["AMIS"]["delta"],
        T=config["AMIS"]["T"],
        target_ess=config["AMIS"]["target_ess_size"],
    script:
        "scripts/sample_parameters.R"


rule sample_parameter_values:
    """
    Samples params["nsamples"] parameter values among weighted values
    generated from the AMIS algorithm.
    Input:
        AMIS CSV output data
    Output:
        CSV data with one column per IU and params["nsamples"] rows.
    """
    input:
        "data/amis_output_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv",
    output:
        "data/sampled_parameters_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv",
    params:
        nsamples=config["nsamples_beta"],
    script:
        "scripts/sample_parameters.py"


rule prepare_mda_file:
    """
    Create input file for forward simulation of infection model.
    Describes start and end simulation year, as well of MDA boundary years.
    Output:
        CSV data

        start_sim_year,end_sim_year,first_mda,last_mda
        2008,2019,2008,2017
    """
    output:
        "data/mda_input.csv",
    script:
        "scripts/prepare_mda_file.py"


rule forward_simulate:
    """
    Simulate infection model for several values of beta parameters, for
    all IUs in a given (first_mda, last_mda, group) subset.
    Input:
        - CSV data with one col per IU and one row per value of
          infection parameter.
        - CSV data describing start and end simulation years, and
          MDA boundary years.
    Output:
        Directory containing model output file for all IUs in subset.
    """
    input:
        "data/sampled_parameters_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv",
        "data/mda_input.csv",
    output:
        directory("data/model_output_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}"),
    script:
        "scripts/resimulate.py"
