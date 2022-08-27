import argparse
from trachoma import Trachoma_Simulation

parser = argparse.ArgumentParser(description="Run trachoma model")
parser.add_argument("--beta-path", required=True)
parser.add_argument("--mda-path", required=True)
parser.add_argument("--pickle-path", required=True)
args = parser.parse_args()

Trachoma_Simulation(
    args.beta_path,
    args.mda_path,
    PrevFilePath="results/PrevFilePath.csv",
    InfectFilePath="results/InfectFilePath.csv",
    SaveOutput=True,
    OutSimFilePath=args.pickle_path,
    logger=None,
)
