import pandas as pd

data = pd.read_csv(snakemake.input[0])

mda_start_years = data.loc[data["start_MDA"] > 0, "start_MDA"]
start_sim_year = mda_start_years.min()

mda_data = pd.DataFrame(
    [[start_sim_year, snakemake.params["end_sim_year"], snakemake.wildcards["FIRST_MDA"], snakemake.wildcards["LAST_MDA"]]],
    columns = ["start_sim_year", "end_sim_year", "first_mda", "last_mda"],
)
slug = "{first_mda}_{last_mda}"
mda_data.to_csv(snakemake.output[0], index=False)
