import argparse
from trachoma import Trachoma_Simulation

parser = argparse.ArgumentParser(description="Run trachoma model")
parser.add_argument("--beta-path", required=True)
parser.add_argument("--mda-path", required=True)
parser.add_argument("--prevalence-path", required=True)
parser.add_argument("--infection-path", required=True)
parser.add_argument("--saved-state", required=True)
args = parser.parse_args()

Trachoma_Simulation(
    args.beta_path,
    args.mda_path,
    PrevFilePath=args.prevalence_path,
    InfectFilePath=args.infection_path,
    SaveOutput=False,
    InSimFilePath=args.saved_state,
    logger=None,
)
