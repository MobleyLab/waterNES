import os, sys
import pickle
import numpy as np
import pathlib
from typing import Dict, List, Optional, Tuple
import warnings
import alchemlyb
sys.path.append('/dfs9/dmobley-lab/swapnilw/waterNES/')
from water_nes.analysis.free_energy_estimate import FreeEnergyEstimate
from water_nes.analysis.nes_free_energy import calculate_nes_free_energy
from alchemlyb.estimators import MBAR
from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.postprocessors.units import get_unit_converter
from alchemlyb.preprocessing import slicing, statistical_inefficiency
from alchemlyb.visualisation import plot_mbar_overlap_matrix
from MDAnalysis import transformations as mda_transformations
from MDAnalysis.analysis import distances as mda_distances
from pymbar.timeseries import ParameterError
from utils import get_water_restraint_force_constant
import matplotlib.pyplot as plt


class Free_Energy:
    def __init__(self, args, stages):
        self.system=args.system
        self.Uscheme=args.Uscheme
        self.Lscheme=args.Lscheme
        self.stages=stages
        self.directory_ligA=args.directory_ligA
        self.directory_ligB=args.directory_ligB
        self.free_energies = {
            "Edge I": FreeEnergyEstimate(value=0, error=0, units="kcal/mol"),
            "Edge J": FreeEnergyEstimate(value=-6.125, error=0.01, units="kcal/mol"),
        }
        if self.directory_ligB == None:
            self.upper_edge(self.directory_ligA)
            #self.lower_edge(self.directory_ligA)
            self.free_energies.update(self.calculate_nes_edges(self.directory_ligA))    # NES Edge
        else:
            self.lower_edge(self.directory_ligB)
        self.write_energies()
        if self.directory_ligA != None:
            self.plot_and_saveA("/dfs9/dmobley-lab/swapnilw/openeye_colab/results")
        if self.directory_ligB != None:
            self.plot_and_saveB("/dfs9/dmobley-lab/swapnilw/openeye_colab/results")

    def calculate_mbar_inner(self,
    stages: List[str],
    cycle_directory: str,
    drop_first: bool,
    file_modifier: Optional[str],
    subsample: bool,
    conservative: bool = True,
    ) -> MBAR:
        u_nk_all_stages = []
        file_name = f"dhdl_{file_modifier}.xvg" if file_modifier is not None else "dhdl.xvg"
        for counter, stage in enumerate(stages):
            with warnings.catch_warnings():
                # Tons of pandas warnings here
                warnings.simplefilter(action="ignore", category=FutureWarning)
                full_u_nk = extract_u_nk(
                    f"{cycle_directory}/stage{stage}/prod/{file_name}", T=298.15
                )
    
                if drop_first:
                    full_u_nk = full_u_nk.drop((0.0, 0.0, 0.0), axis=1)

                u_nk = slicing(full_u_nk, lower=0000)
                if subsample:
                    u_nk_all_stages.append(
                        statistical_inefficiency(
                            u_nk,
                            series=u_nk[u_nk.columns[counter]],
                            conservative=conservative,
                        )
                    )
                else:
                    u_nk_all_stages.append(u_nk)
        #breakpoint()
        return MBAR().fit(alchemlyb.concat(u_nk_all_stages))

    def calculate_mbar(self,
        stages: List[str],
        cycle_directory: str,
        drop_first: bool,
        file_modifier: Optional[str],
        ) -> MBAR:
        """
        TODO: Need to better understand this. MBAR sometimes fails when subsampling,
              and sometimes when _not_ subsampling. The results look reasonable either
              way, so for now I work around this by rerunning the analysis without
              subsampling if needed. This shouldn't affect the current analysis, but it
              would be good to understand the reason for it and possibly fix it!
        """
        try:
            return self.calculate_mbar_inner(
                stages,
                cycle_directory,
                drop_first,
                file_modifier,
                subsample=True,
                conservative=True,
            )
        except ParameterError:
            print("Non-conservative subsampling failed, trying conservative subsampling")
        try:
            return self.calculate_mbar_inner(
                stages,
                cycle_directory,
                drop_first,
                file_modifier,
                subsample=True,
                conservative=False,
            )
        except ParameterError:
            print("Conservative subsampling failed, trying without subsampling")

        return self.calculate_mbar_inner(
            stages, cycle_directory, drop_first, file_modifier, subsample=False
        )

    def calculate_upper_edge_free_energy(self,
    cycle_directory: str, stages: List[str], file_modifier: Optional[str]
) -> Dict[str, FreeEnergyEstimate]:
        output_units = "kcal/mol"
        mbar = self.calculate_mbar(
            stages, cycle_directory, drop_first=False, file_modifier=file_modifier
        )
        do_edge_m = stages == ["1", "4"]
        do_edge_d_endpoint = stages == ["1", "2", "3", "4"]

    # Sanity checks
        if do_edge_m:
            assert mbar.delta_f_.shape == (2, 2)
        elif do_edge_d_endpoint:
            assert mbar.delta_f_.shape == (4, 4)
        else:
            assert mbar.delta_f_.shape[0] > 4 and mbar.delta_f_.shape[1] > 4

        if do_edge_m:
            return {
                "Edge M": FreeEnergyEstimate(
                    value=get_unit_converter(output_units)(mbar.delta_f_).iloc[0][1],
                    error=get_unit_converter(output_units)(mbar.d_delta_f_).iloc[0][1],
                    units=output_units,
                ),
                "upper edge overlap": mbar.overlap_matrix,
            }
        else:
            if do_edge_d_endpoint:
                edge_d = get_unit_converter(output_units)(mbar.delta_f_).iloc[2][3]
                edge_d_error = get_unit_converter(output_units)(mbar.d_delta_f_).iloc[2][3]
            else:
                edge_d = sum(
                    [
                        get_unit_converter(output_units)(mbar.delta_f_).iloc[idx - 1][idx]
                        for idx in range(3, mbar.delta_f_.shape[0])
                    ]
                )
                edge_d_error = sum(
                    [
                        get_unit_converter(output_units)(mbar.d_delta_f_).iloc[idx - 1][idx]
                        for idx in range(3, mbar.d_delta_f_.shape[0])
                    ]
                )

            return {
                "Edge B": FreeEnergyEstimate(
                    value=get_unit_converter(output_units)(mbar.delta_f_).iloc[0][1],
                    error=get_unit_converter(output_units)(mbar.d_delta_f_).iloc[0][1],
                    units=output_units,
                ),
                "Edge C": FreeEnergyEstimate(
                    value=get_unit_converter(output_units)(mbar.delta_f_).iloc[1][2],
                    error=get_unit_converter(output_units)(mbar.d_delta_f_).iloc[1][2],
                    units=output_units,
                ),
                "Edge D": FreeEnergyEstimate(
                    value=edge_d,
                    error=edge_d_error,
                    units=output_units,
                ),
                "upper edge overlap": mbar.overlap_matrix,
            }

    def calculate_lower_edge_free_energy(self,
            cycle_directory: str, stages: List[str], file_modifier: Optional[str]
    ) -> Dict[str, FreeEnergyEstimate]:
        output_units = "kcal/mol"
        restraint_force_constant = get_water_restraint_force_constant(
            cycle_directory=cycle_directory, units="kcal/mol/nm^2"
        )

        rt = 1.9872159e-3 * 298.15  # kcal/mol/K * K
        edge_h = rt * np.log(
            55.5 * 0.6022 * (2 * np.pi * rt / restraint_force_constant) ** (3 / 2)
        )
    # Unit analysis:
    # (2 * np.pi * rt / restraint_force_constant) ** (3 / 2) --> nm^3
    # Molar concentration of water:
    #   55.5 mol / L == 55.5 mol / (10^24 nm^3) == 55.5 * 6.022 * 10^23 / (10^24 nm^3)
    #                == 55.5 * 0.6022 nm^-3
        mbar = self.calculate_mbar(
            stages, cycle_directory, drop_first=True, file_modifier=file_modifier
        )
        do_edge_l = stages == ["7", "5"]
        do_edge_f_endpoint = stages == ["7", "6", "5"]
        #sanity checks
        if do_edge_l:
            assert mbar.delta_f_.shape == (2, 2)
        elif do_edge_f_endpoint:
            assert mbar.delta_f_.shape == (3, 3)
        else:
            assert mbar.delta_f_.shape[0] > 3 and mbar.delta_f_.shape[1] > 3
        if do_edge_l:
            return{
                    "Edge L": FreeEnergyEstimate(
                        value=get_unit_converter(output_units)(mbar.delta_f_).iloc[1][0],
                        error=get_unit_converter(output_units)(mbar.d_delta_f_).iloc[1][0],
                        units=output_units,
                    ),
                    "Edge H": FreeEnergyEstimate(
                        value=edge_h,
                        error=0.0,
                        units=output_units,
                    ),
            "lower edge overlap": mbar.overlap_matrix,
            }
        else:
            edge_f = sum(
                [
                    get_unit_converter(output_units)(mbar.delta_f_).iloc[idx][idx - 1]
                    for idx in range(2, mbar.delta_f_.shape[0])
                ]
            )
            edge_f_error = sum(
                [
                    get_unit_converter(output_units)(mbar.d_delta_f_).iloc[idx][idx - 1]
                    for idx in range(2, mbar.d_delta_f_.shape[0])
                ]
            )

            return {
                "Edge F": FreeEnergyEstimate(
                    value=edge_f,
                    error=edge_f_error,
                    units=output_units,
                ),
                "Edge G": FreeEnergyEstimate(
                    value=get_unit_converter(output_units)(mbar.delta_f_).iloc[1][0],
                    error=get_unit_converter(output_units)(mbar.d_delta_f_).iloc[1][0],
                    units=output_units,
                ),
                "Edge H": FreeEnergyEstimate(
                    value=edge_h,
                    error=0.0,
                    units=output_units,
                ),
                "lower edge overlap": mbar.overlap_matrix,
            }

    def calculate_nes_edges(self, cycle_directory: str) -> Dict[str, FreeEnergyEstimate]:
        edge_and_stages = [
            ("E", 4, 5),
        ]
        num_nes_repeats = 100

        free_energies = {}
        for (
            edge,
            stageA,
            stageB,
        ) in edge_and_stages:
            xvg_files_forward = [
                f"{cycle_directory}/stage{stageA}/NES/run{run}/dhdl.xvg"
                for run in range(1, num_nes_repeats + 1)
            ]
            xvg_files_backward = [
                f"{cycle_directory}/stage{stageB}/NES/run{run}/dhdl.xvg"
                for run in range(1, num_nes_repeats + 1)
            ]
    
            try:
                (
                    free_energies[f"Edge {edge}"],
                    free_energies[f"Work edge {edge}"],
                ) = calculate_nes_free_energy(
                    xvg_files_forward_transition=xvg_files_forward,
                    xvg_files_backward_transition=xvg_files_backward,
                    temperature=298.15,
                    output_units="kcal/mol",
                    bootstrapping_repeats=0,
                )
            except ValueError:
                free_energies[f"Edge {edge}"] = FreeEnergyEstimate(
                    value=np.nan, error=np.nan, units="kcal/mol", bootstrap_error=0
                )
                free_energies[f"Work edge {edge}"] = None
        return free_energies

    def upper_edge(self, cycle_directory):
        self.free_energies.update(
            self.calculate_upper_edge_free_energy(
            cycle_directory=cycle_directory,
            stages=[stage.strip() for stage in self.stages.split(',') if float(stage) < 5],
            file_modifier=f"{self.Uscheme}{self.Lscheme}",
            )
        )
        return

    def lower_edge(self, cycle_directory):
        self.free_energies.update(
            self.calculate_lower_edge_free_energy(
            cycle_directory=cycle_directory,
            stages=[stage.strip() for stage in self.stages.split(',') if float(stage) > 4],
            file_modifier=f"{self.Uscheme}{self.Lscheme}",
            )
        )
        return

    def write_energies(self):
        sum_cycle = 0.0
        error_cycle = 0.0
        for key, value in sorted(self.free_energies.items()):
            if "Edge" in key:
                print(key, value)
                sum_cycle = sum_cycle + value.value
                error_cycle = error_cycle + value.error
        print(f"Edge A   {round(-1 * sum_cycle, 3)} +- {round(error_cycle, 3)} kcal/mol \t from cycle closure")
        with open(f"analysis{self.system}.pickle", "wb") as out_file:
            pickle.dump(self.free_energies, out_file)
        out_file.close()

        return

    def plot_and_saveA(self, out_path):
        fig, axs = plt.subplot_mosaic(
                [["upper", ".", "work", "work"],
                    ["lower", "lower", "work", "work"],
                    ["lower", "lower", "work", "work"]], 
                layout="constrained",
                )
        fig.suptitle(f"NES, Matrix overlap Plots | System {self.system}")

        if "upper edge overlap" in self.free_energies:
            plot_mbar_overlap_matrix(self.free_energies['upper edge overlap'], ax=axs["upper"])

        if "lower edge overlap" in self.free_energies:
            plot_mbar_overlap_matrix(self.free_energies['lower edge overlap'], ax=axs["lower"])

        if "Work edge E" in self.free_energies:
            bins = np.linspace(int(min(np.min(self.free_energies['Work edge E']['forward'] * 0.239006), np.min(self.free_energies['Work edge E']['backward'] * 0.239006))/5) * 5, int(max(np.max(self.free_energies['Work edge E']['forward'] * 0.239006), np.max(self.free_energies['Work edge E']['backward'] * 0.239006))/5) * 5 + 5, 100)
#                    max(self.free_energies['Work edge E']['forward'], self.free_energies['Work edge E']['backward'])/10 + 10, 100)
            axs["work"].hist(self.free_energies['Work edge E']['forward'] * 0.239006, bins, alpha=0.5, label='forward')
            axs["work"].hist(self.free_energies['Work edge E']['backward'] * 0.239006, bins, alpha=0.5, label='backward')
            axs["work"].legend(loc='upper right')
        
        out_file_path = os.path.join(out_path, "_".join(self.system.split('_')[:2]), f"overlap_plots_{self.system.split('_')[0]}_{self.system.split('_')[3]}.png")
        fig.savefig(out_file_path)
        fig.clear()

    def plot_and_saveB(self, out_path):
        fig = plt.figure()
        fig, ax = plt.subplots()
        fig.suptitle(f"Matrix overlap Plot | System {self.system}")
        if "lower edge overlap" in self.free_energies:
            plot_mbar_overlap_matrix(self.free_energies['lower edge overlap'], ax=ax)

        out_file_path = os.path.join(out_path, "_".join(self.system.split('_')[:2]), f"overlap_plots_{self.system.split('_')[1]}_{self.system.split('_')[3]}.png")
        fig.savefig(out_file_path)



