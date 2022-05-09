import pandas as pd
import numpy as np
from numpy.random import default_rng
from pathlib import Path

nsamples = 100
data = pd.read_csv(snakemake.input[0]).set_index("IUCodes")
rng = default_rng()
iu_group = data[
    (data["start_MDA"] == snakemake.wildcards["START_MDA"]) &&
    (data["last_MDA"] == snakemake.wildcards["LAST_MDA"])  &&
    (data["group"] == snakemake.wildcards["GROUP"])
]
prevs = iu_group.apply(
    lambda s: rng.normal(s["Logit"], s["Sds"], nsamples),
    axis=1,
    result_type="expand"
).apply(lambda x: np.exp(x) / (1 + np.exp(x)))
prevs.to_csv(snakemake.output[0], index=False)

