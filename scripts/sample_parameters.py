import pandas as pd
import numpy as np
from numpy.random import default_rng

data = pd.read_csv(snakemake.input[0])
rng = default_rng()
iu_weights = data.drop(["seeds", "beta", "sim_prev"], axis=1)
nsamples = snakemake.params["nsamples"]
for col in iu_weights:
    df = data[["seeds", "beta"]].sample(
        nsamples,
        weights=iu_weights[col],
        replace=False
    ).to_csv(snakemake.output[0], index=False)

