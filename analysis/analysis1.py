import argparse
from utils import * 
from reading_file import Reading_analysis
from processing_data import Processing_Data
from calculate_free_energy import Free_Energy
from plot_distances import plot_dists


def input_data():
    parser = argparse.ArgumentParser(description="Analyze ABFE/RBFE for system")
    parser.add_argument("--system", type=str, help="The name of the system")
    parser.add_argument("--fe", type=str, help="Type of calculation: abfe or rbfe")
    parser.add_argument(
        "--directory_ligA",
        type=str,
        help=(
            "The top directory of the cycle to analyze for Ligand A, should contain folders "
            "named stageX where X is one of the valid stage numbers"
        ),
    )
    parser.add_argument(
        "--directory_ligB",
        type=str,
        help=(
            "The top directory of the cycle to analyze for Ligand B, should contain folders "
            "named stageX where X is one of the valid stage numbers. Only required for RBFE."
        ),
        required=False
    )
    parser.add_argument("--dist", action="store_true", help="Do attachment site- buried water distance calculations", required=False, default=False)
    parser.add_argument(
        "--postProcessTrapped",
        action="store_true",
        help="Postprocess the trapped water and solvent restraint in stage 1",
        required=False
    )
    parser.add_argument(
        "--postProcessPocket",
        action="store_true",
        help="Postprocess the binding pocket restraints",
        required=False
    )
    parser.add_argument("--Uscheme", 
            type=str, 
            help="Scheme for upper edge: 1 for 1-4, 2 for 1-2-3-4, 3 for 1-2-3-3.1-3.2...-4", 
            required=False, default="2")
    parser.add_argument("--Lscheme", 
            type=str,
            help="Scheme for lower edge: 1 for 5-7, 2 for 5-6-6.1-6.2-....-7", 
            required=False, default="1")
    parser.add_argument("--read", action="store_true",
            help="Just read the pickle file that stores the analysis",
            required=False, default=False)
    parser.add_argument("--ligand", action="store_true",
            help="Name of the ligand residue. Default=LIG",
            required=False, default="LIG")
    args = parser.parse_args()
    return args
    

def calculate_free_energies(args, stages):
    free_energies=Free_Energy(args, stages)
    return

def reading_data(args):
    reading_analysis = Reading_analysis(args)
    sanity_check = reading_analysis.sanity_check()
    assert(sanity_check)
    return reading_analysis

def processing_data(reading_analysis):
    processed_data = Processing_Data(reading_analysis)
    #processed_data.engine()
    return processed_data

def main():
    args = input_data()
    reading_analysis = reading_data(args)
    processing_data(reading_analysis)
    calculate_free_energies(args, reading_analysis.stages)
    #plot_dists(reading_analysis.metadata)
    return

if __name__ == "__main__":
    main()


