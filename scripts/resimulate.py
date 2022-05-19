import pandas as pd
import os
from pathlib import Path
from trachoma import Trachoma_Simulation

slug = "_".join(snakemake.wildcards)
pin = Path(f"data/model_input_{slug}")
pin.mkdir(parents=True, exist_ok=True)
pout = Path(f"data/model_output_{slug}")
pout.mkdir(parents=True, exist_ok=True)
    
sampled_params = pd.read_csv(snakemake.input[0])
for iucode in sampled_params.columns:
    beta_input = pin / f"beta_values_{iucode}.csv"
    # Using dataframe index as seed value
    sampled_params[iucode].to_csv(beta_input, header=True, index=True)
    # Have to convert paths to strings because simulation methods
    # iterate on filename strings at some point
    Trachoma_Simulation(
        str(beta_input),
        snakemake.input[1],
        str(pout / f"prevalence_{iucode}.csv"),
        str(pout / f"infection_{iucode}.csv"),
        SaveOutput=True,
        OutSimFilePath=str(pout / f"output_vals_{iucode}.p"),
        InSimFilePath=None
    )
