configfile: "config.yaml"

def aggregate_input(wildcards):
    from pandas import read_csv
    data = read_csv(config["data"])
    return [f"results/infection_{iucode}.csv" for iucode in data["IUCodes"]]

localrules: all, group_ius, make_mda_file, prepare_mda_file

rule all:
    input:
        aggregate_input,


checkpoint group_ius:
    """
    Assign each iu a prevalence level according to mean prevalence value.
    For example
    level 1: 0 < prev <= 10%
    level 2: 10 < prev <= 20%
    ...
    input:
        CSV data describing one IU per row, and columns 'Logit' and 'Sds'
        describing logit and standard deviation values for each IU.
    output:
        Input CSV data augmented with extra column 'level'.
    """
    input:
        config["data"],
    output:
        "results/FinalDataLevel.csv",
    script:
        "scripts/group_ius.py"


rule make_prevalence_maps:
    """
    Sample params['nsamples'] prevalence values for each IU in
    a (first_mda, last_mda, group) group.
    input:
        CSV data describing one IU per row:

        Logit,Sds,start_MDA,last_MDA,level
        -1.814302776,0.243590266,2008,2017,level_2
        ...

    output:
        CSV data with one row per IU in (first_mda, last_mda, group) group
        and params["nsamples"] columns.
    """
    input:
        "results/FinalDataLevel.csv",
    output:
        "results/prev_map_{FIRST_MDA}_{LAST_MDA}_level_{LEVEL}.csv",
    params:
        nsamples=config["nsamples_prevalence_map"],
    group: "amis"
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
        config["data"],
    output:
        "results/mda_input_{FIRST_MDA}_{LAST_MDA}.csv",
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
        "results/mda_input_{FIRST_MDA}_{LAST_MDA}.csv",
        "results/prev_map_{FIRST_MDA}_{LAST_MDA}_level_{LEVEL}.csv",
    output:
        "results/amis_output_{FIRST_MDA}_{LAST_MDA}_level_{LEVEL}.csv",
    params:
        delta=config["AMIS"]["smooth_param"],
        nsamples=config["AMIS"]["nsamples"],
        mixture_samples=config["AMIS"]["mixture_samples"],
        df=config["AMIS"]["df"],
        target_ess=config["AMIS"]["target_ess_size"],
        log=config["AMIS"]["log_scale"],
        max_iters=config["AMIS"]["max_iterations"],
    group: "amis"
    script:
        "scripts/sample_parameters.R"


def get_subset_from_IU_code(wildcards):
    from pandas import read_csv
    checkpoint_output = checkpoints.group_ius.get(**wildcards).output[0]
    iucode = wildcards["IUCODE"]
    row = read_csv(checkpoint_output).set_index("IUCodes").loc[iucode]
    return (row["start_MDA"], row["last_MDA"], row["level"])


def get_input_amis_file(wildcards):
    start_MDA, last_MDA, level = get_subset_from_IU_code(wildcards)
    return f"results/amis_output_{start_MDA}_{last_MDA}_{level}.csv"


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
        get_input_amis_file,
    output:
        "results/sampled_parameters_{IUCODE}.csv",
    params:
        nsamples=config["nsamples_beta"],
    group: "resimulate"
    script:
        "scripts/sample_parameters.py"


def get_input_mda_file(wildcards):
    start_MDA, last_MDA, _ = get_subset_from_IU_code(wildcards)
    return f"results/mda_input_{start_MDA}_{last_MDA}.csv"

rule resimulate_history:
    input:
        "results/sampled_parameters_{IUCODE}.csv",
        get_input_mda_file,
    output:
        "results/output_state_{IUCODE}.p"
    group: "resimulate"
    shell:
        """
        python scripts/resimulate_history.py \
            --beta-path {input[0]} \
            --mda-path {input[1]} \
            --pickle-path {output[0]}
        """

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
        "results/mda_input.csv",
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
        sampled_parameters="results/sampled_parameters_{IUCODE}.csv",
        saved_state="results/output_state_{IUCODE}.p",
        mda_input="results/mda_input.csv",
    output:
        infection="results/infection_{IUCODE}.csv",
        prevalence="results/prevalence_{IUCODE}.csv",
    group: "resimulate"
    run:
        from trachoma import Trachoma_Simulation

        Trachoma_Simulation(
            input.sampled_parameters,
            input.mda_input,
            PrevFilePath=output.prevalence,
            InfectFilePath=output.infection,
            SaveOutput=False,
            InSimFilePath=input.saved_state,
            logger=None,
        )
