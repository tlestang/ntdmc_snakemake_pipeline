from trachoma import Trachoma_Simulation

Trachoma_Simulation(
  snakemake.input["sampled_parameters"],
  snakemake.input["mda_input"],
  PrevFilePath=snakemake.output["prevalence"],
  InfectFilePath=snakemake.output["infection"],
  SaveOutput=False,
  InSimFilePath=snakemake.input["saved_state"],
  logger=None,
  num_cores=snakemake.threads,
)
