import pandas as pd
import numpy as np
from numpy.random import default_rng

data = pd.read_csv(snakemake.input[0])
rng = default_rng()

cols_to_drop = list(data.filter(regex = "prev"))
cols_to_drop += ["seeds", "beta"]
iu_weights = data.drop(cols_to_drop, axis=1)

nsamples = snakemake.params["nsamples"]
for col in iu_weights:
    df = data[["seeds", "beta"]].sample(
        nsamples,
        weights=iu_weights[col],
        replace=False
    ).to_csv(snakemake.output[0], index=False)

