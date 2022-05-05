import pandas as pd
import os
from trachoma import Trachoma_Simulation

pin = Path("data/model_input")
pin.mkdir(parents=True, exist_ok=True)
pout = Path("data/model_output")
pout.mkdir(parents=True, exist_ok=True)
    
sampled_params = pd.read_csv("snakemake.input[0]")
for iucode in sampled_params.columns:
    beta_input = pin / f"beta_values_{iucode}.csv"
    # Using dataframe index as seed value
    sampled_params[iucode].to_csv(beta_input, header=True, index=True)
    with open(pout / f"prevalence_{iucode}.csv") as f:
        f.write("prevalence output here")
    with open(pout / f"infection_{iucode}.csv") as f:
        f.write("infection output here")
    # Trachoma_Simulation(
    #     beta_input,
    #     snakemake.input[1],
    #     pout / f"prevalence_{iucode}.csv",
    #     pout / f"infection_{iucode}.csv",
    #     SaveOutput=True,
    #     OutSimFilePath=pout / f"output_vals_{iucode}.p"
    #     InSimFilePath=None
    # )
