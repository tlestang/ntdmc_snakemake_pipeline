import pandas as pd
import numpy as np
from numpy.random import default_rng

data = pd.read_csv(snakemake.input[0])
rng = default_rng()
iu_weights = data.drop(["seeds", "beta", "sim_prev"], axis=1)
beta_values = data["beta"]
nsamples = snakemake.params["nsamples"]
sampled_params = iu_weights.apply(
    lambda x: np.random.choice(beta_values, (nsamples,), replace=False, p=x),
    axis=0,
    result_type="expand"
)
sampled_params.to_csv(snakemake.output[0], index=False)

