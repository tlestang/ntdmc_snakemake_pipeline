import pandas as pd
import numpy as np
from numpy.random import default_rng
from pathlib import Path

nsamples = 100
data = pd.read_csv(snakemake.input[0]).set_index("IUCodes")
rng = default_rng()
for name, group in data.groupby(["start_MDA", "last_MDA", "group"]):
    prevs = group.apply(
        lambda s: rng.normal(s["Logit"], s["Sds"], nsamples),
        axis=1,
        result_type="expand"
    ).apply(lambda x: np.exp(x) / (1 + np.exp(x)))
    slug = f"{name[0]}_{name[1]}_{name[2]}"
    p = Path("data/prevalence_maps")
    p.mkdir(parents=True, exist_ok=True)
    prevs.to_csv(p / f"prev_map_{slug}.csv", index=False)
