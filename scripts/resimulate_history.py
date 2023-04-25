from trachoma import Trachoma_Simulation

Trachoma_Simulation(
    snakemake.input[0],
    snakemake.input[1],
    PrevFilePath="results/PrevFilePath.csv",
    InfectFilePath="results/InfectFilePath.csv",
    SaveOutput=True,
    OutSimFilePath=snakemake.output[0],
    logger=None,
    num_cores=snakemake.threads,
)
