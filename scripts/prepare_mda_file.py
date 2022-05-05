import os
import pandas as pd

mda_data = pd.DataFrame(
    [[2007, 2019, 2008, 2017]],
    columns = ["start_sim_year", "end_sim_year", "first_mda", "last_mda"],
)
mda_data_input = os.path.join(snakemake.output[0])
mda_data.to_csv(mda_data_input, index=False)
