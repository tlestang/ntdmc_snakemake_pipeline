import pandas as pd
import numpy as np
from numpy.random import default_rng
from pathlib import Path

nsamples = snakemake.params["nsamples"]
data = pd.read_csv(snakemake.input[0]).set_index("IUCodes")
rng = default_rng()
iu_group = data[
    (data["start_MDA"] == int(snakemake.wildcards["FIRST_MDA"])) &
    (data["last_MDA"] == int(snakemake.wildcards["LAST_MDA"]))  &
    (data["level"] == f"level_{snakemake.wildcards['LEVEL']}")
]
prevs = iu_group.apply(
    lambda s: rng.normal(s["Logit"], s["Sds"], nsamples),
    axis=1,
    result_type="expand"
).apply(lambda x: np.exp(x) / (1 + np.exp(x)))
prevs.to_csv(snakemake.output[0], index=False)

