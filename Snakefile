configfile: "config.yaml"

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.make_prevalence_maps.get(**wildcards).output[0]
    wcards = glob_wildcards(os.path.join(checkpoint_output, "prev_map_{slug}.csv"))
    return expand(
        os.path.join("data", "amis_output_{slug}.csv"), slug=wcards.slug
    )

rule all:
    input: aggregate_input
    shell:
        "echo INPUT IS {input}"

rule make_group_scenario_pairs:
    input:
        "data/FinalDataTest.csv",
    output:
        "data/FinalDataGroup.csv"
    script:
        "scripts/group_ius.py"

checkpoint make_prevalence_maps:
    input:
        "data/FinalDataGroup.csv"
    output:
        directory("data/prevalence_maps")
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

def prevalence_map_input(wildcards):
    checkpoint_output = checkpoints.make_prevalence_maps.get(**wildcards).output[0]
    return os.path.join(checkpoint_output, "prev_map_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv")


rule estimate_parameter_weights:
    input:
        "data/mda_input_{FIRST_MDA}_{LAST_MDA}.csv",
        prevalence_map_input
    output:
        "data/amis_output_{FIRST_MDA}_{LAST_MDA}_group_{GROUP}.csv"
    params:
        nsamples=100,
        delta=5,
        T=2,
        target_ess=250
    script:
        "scripts/sample_parameters.R"
