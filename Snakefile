configfile: "config.yaml"

def aggregate_input(wildcards):
    from pandas import read_csv
    data = read_csv(config["data"])
    return [f"results/infection_{iucode}.csv" for iucode in data["IUCodes"]]

localrules: all, group_ius


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

rule the_last_rule:
    input:
        "results/prev_map_{FIRST_MDA}_{LAST_MDA}_level_{LEVEL}.csv"
    output:
        "results/output_{FIRST_MDA}_{LAST_MDA}_level_{LEVEL}.csv"
    group: "amis"
    script:
        "scripts/dummy.py"
