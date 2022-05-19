configfile: "config.yaml"

rule forward_simulate:
    input:
        "data/sampled_parameters_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv",
        "data/mda_input.csv"
    output:
        directory("data/model_output_{FIRST_MDA}_{LAST_MDA}_{GROUP}")
    script:
        "scripts/resimulate.py"


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
    input:
        "data/FinalDataGroup.csv"
    output:
        "data/prev_map_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv"
    script:
        "scripts/compute_prevalence_maps.py"

rule make_mda_file:
    input:
        "data/FinalDataTest.csv",
    output:
        "data/mda_input_{FIRST_MDA}_{LAST_MDA}.csv"
    params:
        END_SIM_YEAR=2019
    script:
        "scripts/make_mda_files.py"

rule estimate_parameter_weights:
    input:
        "data/mda_input_{FIRST_MDA}_{LAST_MDA}.csv",
        "data/prev_map_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv"
    output:
        "data/amis_output_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv"
    params:
        nsamples=100,
        delta=5,
        T=2,
        target_ess=250
    script:
        "scripts/sample_parameters.R"

rule sample_parameter_values:
    input:
        "data/amis_output_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv"
    output:
        "data/sampled_parameters_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv"
    params:
        nsamples=10
    script:
        "scripts/sample_parameters.py"

rule prepare_mda_file:
    output:
        "data/mda_input.csv"
    script:
        "scripts/prepare_mda_file.py"
