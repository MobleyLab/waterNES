import os, sys
import pickle
import pathlib
import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations as mda_transformations
from MDAnalysis.analysis import distances as mda_distances
from utils import calculate_distances, get_pocket_selection_string, get_pocket_restraint_force_constant, simple_restraint, get_water_restraint_force_constant, get_solvent_restraint_c12, vdw_repulsive_energy
from plot_distances import plot_dists


class Reading_analysis:
    def __init__(self, args):
        print("reading the information...")
        #self.fe_option=args.fe_option
        #self.cycle_directory=args.directory
        self.system=args.system
        self.fe = args.fe
        self.cycle_directory = None
        self.directory_ligA = args.directory_ligA
        self.directory_ligB = args.directory_ligB
        if args.directory_ligB == None:
            self.cycle_directory = args.directory_ligA
        elif args.directory_ligA == None:
            self.cycle_directory = args.directory_ligB
        self.Uscheme = args.Uscheme
        self.Lscheme = args.Lscheme
        self.read=args.read
        self.postprocess_trapped=args.postProcessTrapped
        self.postprocess_pocket=args.postProcessPocket
        self.ligand = args.ligand
        self.stages = self.determine_stages()
        self.analysis = {}
        self.metadata, self.energies = self.engine()

    def add_workflow(self, universe):
        protein_and_ligand = universe.select_atoms("protein or resname %s" % self.ligand)
        water_and_virtual_site = universe.select_atoms("resname MOL or resname ATT")

        center_selection = universe.select_atoms("resname %s" % self.ligand)
        if len(center_selection) == 0:
            center_selection = universe.select_atoms("resname MOL")

        assert len(center_selection) > 0

        workflow = [
            mda_transformations.unwrap(universe.atoms),
            mda_transformations.center_in_box(center_selection, center="mass"),
            mda_transformations.wrap(protein_and_ligand, compound="fragments"),
            mda_transformations.wrap(water_and_virtual_site, compound="atoms"),
        ]
        universe.trajectory.add_transformations(*workflow)
        return universe

    def load_universe(self, topology, trajectory, transfer_to_memory=True, step=1, workflow=True):
        universe = mda.Universe(topology, trajectory)
        if transfer_to_memory:
            universe.transfer_to_memory(step=step)
        if workflow:
            return self.add_workflow(universe)
        return universe


    def distance_calcs(self, all_stages, pocket_selection_string):
        distances = {}
        if pathlib.Path(f"analysis{self.system}_dist.pickle").is_file():
            with open(f"analysis{self.system}_dist.pickle", "rb") as in_file:
                self.analysis["distances"] = pickle.load(in_file)
                distances = self.analysis["distances"]
            in_file.close()
        for stage in all_stages:
            if f"stage{stage}" in distances.keys():
                continue
            else:
                trajectory = self.load_universe(
                        topology=f"{self.cycle_directory}/stage{stage}/prod/topol.tpr",
                        trajectory=f"{self.cycle_directory}/stage{stage}/prod/traj_comp.xtc",
                    )
                distances[f"stage{stage}"] = calculate_distances(
                        universe=trajectory, pocket_selection_string=pocket_selection_string
                    )
            self.analysis["distances"] = distances
        with open(f"analysis{self.system}_dist.pickle", "wb") as out_file:
            pickle.dump(self.analysis["distances"], out_file)
        out_file.close()
        return distances

    def trapped_calcs(self):
        trapped_restraint_energy = {}
        solvent_restraint_energy = {}
        if self.postprocess_trapped:
            restraint_force_constant = get_water_restraint_force_constant(
                cycle_directory=self.cycle_directory, units="kJ/mol/Å^2"
            )
            trapped_restraint_energy["stage1"] = simple_restraint(
                force_constant=restraint_force_constant,
                distance=self.analysis["distances"]["stage1"]["closest water"],
            )
            restraint_c12 = get_solvent_restraint_c12(
                cycle_directory=self.cycle_directory, units="kJ/mol*Å^12"
            )
            solvent_restraint_energy_old = vdw_repulsive_energy(
                c12=restraint_c12,
                distance=self.analysis["distances"]["stage1"]["solvent water 1st"],
            )
            solvent_restraint_energy_old += vdw_repulsive_energy(
                c12=restraint_c12,
                distance=self.analysis["distances"]["stage1"]["solvent water 2nd"],
            )
            solvent_restraint_energy_new = vdw_repulsive_energy(
                c12=restraint_c12,
                distance=self.analysis["distances"]["stage1"]["closest solvent"]
            )
            solvent_restraint_energy_new += vdw_repulsive_energy(
                c12=restraint_c12,
                distance=self.analysis["distances"]["stage1"]["second closest solvent"]
            )
            solvent_restraint_energy_new += vdw_repulsive_energy(
                c12=restraint_c12,
                distance=self.analysis["distances"]["stage1"]["third closest solvent"]
            ) 
            solvent_restraint_energy_new += vdw_repulsive_energy(
                c12=restraint_c12,
                distance=self.analysis["distances"]["stage1"]["fourth closest solvent"]
            ) 
            solvent_restraint_energy_new += vdw_repulsive_energy(
                c12=restraint_c12,
                distance=self.analysis["distances"]["stage1"]["fifth closest solvent"]
            )
            solvent_restraint_energy["stage1"] = np.array(
                [solvent_restraint_energy_old, solvent_restraint_energy_new]
            )
        return trapped_restraint_energy, solvent_restraint_energy 
            
    def pocket_restraints_calcs(self, all_non_restraint_stages, pocket_selection_string, reference):
        position_restraint_energy = {}
        if self.postprocess_pocket:
            if os.path.isfile(f"analysis{self.system}_ppocket.pickle"):
                with open(f"analysis{self.system}_ppocket.pickle", "rb") as in_file:
                    position_restraint_energy = pickle.load(in_file)
                in_file.close()
            for stage in all_non_restraint_stages:
                if f"stage{stage}" in position_restraint_energy.keys():
                    continue
                else:
                    trajectory = self.load_universe(
                        topology=f"{self.cycle_directory}/stage{stage}/prod/topol.tpr",
                        trajectory=f"{self.cycle_directory}/stage{stage}/prod/traj_comp.xtc",
                        transfer_to_memory=False,
                    )
                # Align frames, and calculate positional restraint
                    position_restraint_energy[f"stage{stage}"] = np.zeros(
                        len(trajectory.trajectory)
                    )
                    force_constant = get_pocket_restraint_force_constant(
                        cycle_directory=self.cycle_directory, units="kJ/mol/Å^2"
                    )
                    for (
                        frame_idx,
                        _,
                    ) in enumerate(trajectory.trajectory):
                        mda.analysis.align.alignto(
                            trajectory.select_atoms(pocket_selection_string),
                            reference.select_atoms(pocket_selection_string),
                        )
                        distances = np.linalg.norm(
                            trajectory.select_atoms(pocket_selection_string).positions
                            - reference.select_atoms(pocket_selection_string).positions,
                            axis=1,
                        )
                        position_restraint_energy[f"stage{stage}"][
                            frame_idx
                        ] = simple_restraint(
                            force_constant=force_constant,
                            distance=distances,
                        ).sum()
            with open(f"analysis{self.system}_ppocket.pickle", "wb") as out_file:
                pickle.dump(position_restraint_energy, out_file)
            out_file.close()
        return position_restraint_energy

    def make_analysis(self):
        trapped_restraint_energy = None
        solvent_restraint_energy = None
        reference = None
        if self.directory_ligA != None:
            reference = self.load_universe(
                topology=f"{self.cycle_directory}/stage1/min/topol.tpr",
                trajectory=f"{self.cycle_directory}/../minimized.gro",
            )
        else:
            reference = self.load_universe(
                topology=f"{self.cycle_directory}/stage5/min/topol.tpr",
                trajectory=f"{self.cycle_directory}/../minimized.gro",
            )
        pocket_selection_string = get_pocket_selection_string(reference)
        all_stages= [
                stage.strip()
                for stage in self.stages.split(',')
        ]
        all_non_restraint_stages=[
                stage
                for stage in all_stages if stage != "4" and stage != "5"
        ]
        distances = self.distance_calcs(all_stages, pocket_selection_string)
        if self.directory_ligA != None:
            plot_dists(distances, self.system, "A")
        else:
            plot_dists(distances, self.system, "B")
        if self.directory_ligA != None:
            trapped_restraint_energy, solvent_restraint_energy = self.trapped_calcs()
        pocket_restraints = self.pocket_restraints_calcs(all_non_restraint_stages, pocket_selection_string, reference)
        return distances, trapped_restraint_energy, solvent_restraint_energy, pocket_restraints


    def engine(self):
        metadata = {}
        energies = {}
        if self.read:
            with open(f"analysis{self.system}.pickle", "rb") as in_file:
                self.analysis = pickle.load(in_file)
            in_file.close()
        else:
            metadata["distances"], energies["trapped_restraint_energy"], energies["solvent_restraint_energy"], energies["pocket_restraints"] = self.make_analysis()
#        return self.stages, self.trapped_restraint_energy, self.solvent_restraint_energy, self.position_restraint_energy, self.analysis
        #return self.stages, trapped_restraint_energy, solvent_restraint_energy, pocket_restraints, self.analysis
        return metadata, energies

    def sanity_check(self):
        status = True
        for stage in self.stages.split(','):
            if os.path.isfile(f"{self.cycle_directory}/stage{stage.strip()}/min/topol.tpr") == False:
                status = False
                print(f"missing file: {self.cycle_directory}/stage{stage.strip()}/min/topol.tpr")
            elif os.path.isfile(f"{self.cycle_directory}/stage{stage.strip()}/prod/traj_comp.xtc") == False:
                status = False
                print(f"missing file: {self.cycle_directory}/stage{stage.strip()}/prod/traj_comp.xtc")

        return status

    def determine_stages(self):
        schemes_to_Ustages = {"1": "1, 4", "2": "1, 2, 3, 4", "3": "1, 2, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 4"}
        schemes_to_Lstages = {"1": "7, 5", "2": "7, 6, 5", "3": "7, 6, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1, 5"}
        Ustages = ""
        Lstages = ""
        stages = ""
        if self.Uscheme != None:
            Ustages = schemes_to_Ustages[self.Uscheme]
            stages = stages + Ustages 
        if self.Lscheme != None:
            Lstages = schemes_to_Lstages[self.Lscheme]
            if stages == "":
                stages = Lstages
            else:
                stages = stages + "," + Lstages

        #return fe_option_to_stages[self.fe_option]
        return stages


