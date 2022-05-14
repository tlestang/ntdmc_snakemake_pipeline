import pandas as pd
from math import ceil
import numpy as np
from numpy.random import default_rng

def estimate_mean_prevalence(loc, scale, nsamples):
    rng = default_rng()
    samples = rng.normal(loc, scale, size=nsamples)
    sample_prevs = np.exp(samples) / (1 + np.exp(samples))
    return np.mean(sample_prevs) * 100


def comp_prevalence_group(iu_data, nsamples=3000):
    prev = estimate_mean_prevalence(iu_data["Logit"], iu_data["Sds"], nsamples)
    if prev >= 60:
        return 7
    return f"group_{ceil(prev / 10)}"

data = pd.read_csv(snakemake.input[0])
data["group"] = data.apply(comp_prevalence_group, axis=1, nsamples=3000)
data.to_csv(snakemake.output[0], index=False)
