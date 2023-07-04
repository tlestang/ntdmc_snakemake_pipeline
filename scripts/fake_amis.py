import csv
import pandas as pd


FIRST_MDA = snakemake.wildcards["FIRST_MDA"]
LAST_MDA = snakemake.wildcards["LAST_MDA"]
LEVEL = snakemake.wildcards["LEVEL"]
NSAMPLES = snakemake.params["nsamples"]

df = pd.read_csv(snakemake.input[2]).set_index('IUCodes')
IU_subset = df.groupby(
    ['start_MDA', 'last_MDA', 'level']
).get_group((int(FIRST_MDA), int(LAST_MDA), f'level_{LEVEL}'))


fname = f"results/amis_output_{FIRST_MDA}_{LAST_MDA}_level_{LEVEL}.csv"
with open(fname, 'w', newline='') as f:
    csvwriter = csv.writer(f, delimiter=',')
    headerrow = (
        ["seeds", "beta", "prev"] +
        [iucode for iucode in IU_subset.index]
    )
    csvwriter.writerow(headerrow)
    for i in range(NSAMPLES):
        csvwriter.writerow(
            [1, 0.12, 0.4] + [1.] * len(IU_subset.index)
        )
