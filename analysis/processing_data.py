import pandas as pd
import os, sys
import pathlib
import numpy as np
np.set_printoptions(suppress=True)
from typing import Optional, Tuple, Dict, List
from check_water_positions import Water_Positions_Analysis
from dhdl_legends import Custom_Legends

class Processing_Data:
    def __init__(self, reading_analysis):
        print("processing the data ...")
        self.reading_analysis = reading_analysis
        self.status = self.engine()


    def water_position_analysis(self, current_stage):
        mask = np.full((601), False)
        if f"stage{current_stage}" in self.reading_analysis.metadata['distances']:

            print(f"checking water positions for stage {current_stage}")
            trapped_water = self.reading_analysis.metadata['distances'][f"stage{current_stage}"]['trapped water']
            closest_water = self.reading_analysis.metadata['distances'][f"stage{current_stage}"]['closest water']
            closest_solvent = self.reading_analysis.metadata['distances'][f"stage{current_stage}"]['closest solvent']
            first_solvent = self.reading_analysis.metadata['distances'][f"stage{current_stage}"]['solvent water 1st']
            eq_dist = None
            if 'stage1' in self.reading_analysis.metadata['distances']:
                eq_dist = np.average(self.reading_analysis.metadata['distances']['stage1']['trapped water'][0:5])
            else:
                eq_dist = np.average(self.reading_analysis.metadata['distances']['stage5']['trapped water'][0:10])
#            eq_dist = 1.4
            print(eq_dist)
            if int(float(current_stage)) <= 4:
                mask = np.logical_and(closest_water <= eq_dist + 1.5, closest_solvent >= closest_water + 1)
            elif int(float(current_stage)) >= 5:
                mask = np.logical_and(1, first_solvent >= eq_dist)
            print(f"{round(mask.sum()/len(closest_water) * 100, 3)}% data used")
#            breakpoint()
        else:
            print(f"water positions for stage {current_stage} is not calculated; taking all data from this stage")
            mask = np.full((601), True)

        return mask 

    def read_xvg(self, current_stage, xvg_file_name: str) -> Tuple[Dict[str, np.ndarray], List[str]]:
        contents = {}
        header = []
        with open(xvg_file_name, "r") as xvg_file:
            for line in xvg_file:
                if (line.startswith("#") or line.startswith("@")):
                    if "legend" not in line:
                        header.append(line.strip())
                        continue
                    else:
                        continue
                entries = line.split()
                time = entries[0]
                #assert time not in contents
                contents[time] = np.array([float(entry) for entry in entries[1:]])
        legends=Custom_Legends()
        if int(float(current_stage)) < 5:
            if self.reading_analysis.Uscheme == "1":
                legend = legends.legend3
            elif self.reading_analysis.Uscheme == "2":
                legend = legends.legend2
            elif self.reading_analysis.Uscheme == "3":
                legend = legends.legend1
        elif int(float(current_stage)) > 4:
            if self.reading_analysis.Lscheme == "1":
                legend = legends.legend4
            elif self.reading_analysis.Lscheme == "2":
                legend = legends.legend2
            elif self.reading_analysis.Lscheme == "3":
                legend = legends.legend1
        return contents, header, legend

    def write_xvg(self, current_stage) -> None:
        #mod_file = pathlib.Path(f"{self.reading_analysis.cycle_directory}/stage{stage}/prod/dhdl_{self.reading_analysis.Uscheme}{self.reading_analysis.Lscheme}.xvg")
        input_xvg_file=pathlib.Path(f"{self.reading_analysis.cycle_directory}/stage{current_stage}/prod/dhdl.xvg")
        output_xvg_file=pathlib.Path(f"{self.reading_analysis.cycle_directory}/stage{current_stage}/prod/dhdl_{self.reading_analysis.Uscheme}{self.reading_analysis.Lscheme}.xvg")
        input_xvg_contents, header, legend = self.read_xvg(current_stage, xvg_file_name=input_xvg_file)
        num_entries = len(input_xvg_contents)
        if self.reading_analysis.energies["trapped_restraint_energy"] != None and f"stage{current_stage}" in self.reading_analysis.energies["trapped_restraint_energy"]:
            assert len(self.reading_analysis.energies["trapped_restraint_energy"][f"stage{current_stage}"]) == num_entries
            assert current_stage == "1"
        if self.reading_analysis.energies["solvent_restraint_energy"] != None and f"stage{current_stage}" in self.reading_analysis.energies["solvent_restraint_energy"]:
            assert self.reading_analysis.energies["solvent_restraint_energy"][f"stage{current_stage}"].shape == (2, num_entries)
            assert current_stage == "1"
        if f"stage{current_stage}" in self.reading_analysis.energies["pocket_restraints"]:
            assert len(self.reading_analysis.energies["pocket_restraints"][f"stage{current_stage}"]) == num_entries
            assert current_stage not in ["4", "5"]

        reference_to_column = {'1': '000',
                                '2': '010',
                                '3': '110',
                                '4': '111',
                                '5': '111',
                                '6': '110',
                                '7': '010',
                                '3.1': '11p01',
                                '3.2': '11p02',
                                '3.3': '11p04',
                                '3.4': '11p08',
                                '3.5': '11p16',
                                '3.6': '11p32',
                                '3.7': '11p64',
                                '6.1': '11p01',
                                '6.2': '11p02',
                                '6.3': '11p04',
                                '6.4': '11p08',
                                '6.5': '11p16',
                                '6.6': '11p32',
                                '6.7': '11p64'
                                }

        df_xvg_file = pd.DataFrame.from_dict(input_xvg_contents, orient='index')
        df_xvg_file = df_xvg_file.reset_index()
        df_xvg_file.columns = ['time', 'PE', 'vdw_l', 'bonded_l', 'restraint_l', '000', '010', '110', '11p01', '11p02', '11p04', '11p08', '11p16', '11p32', '11p64', '111', 'pV']
        if f"stage{current_stage}" in self.reading_analysis.energies["pocket_restraints"]:
            df_xvg_file.drop('restraint_l', axis=1, inplace=True)
            df_pocket_restraint_energy = pd.DataFrame(self.reading_analysis.energies["pocket_restraints"][f"stage{current_stage}"], columns=['restraint_l'])
            df_xvg_file = pd.concat([df_xvg_file, df_pocket_restraint_energy], axis=1)
            if current_stage == '1':
                df_xvg_file.drop('bonded_l', axis=1, inplace=True)
                df_trapped_water_energy = pd.DataFrame(self.reading_analysis.energies["trapped_restraint_energy"][f"stage{current_stage}"], columns=['bonded_l'])
                df_xvg_file = pd.concat([df_xvg_file, df_trapped_water_energy], axis=1)
                df_xvg_file.drop('vdw_l', axis=1, inplace=True)
                df_solvent_restraint_energy = pd.DataFrame(np.transpose(self.reading_analysis.energies["solvent_restraint_energy"][f"stage{current_stage}"]), columns=['solvent_restraint_energy_old', 'vdw_l'])
                df_xvg_file = pd.concat([df_xvg_file, df_solvent_restraint_energy], axis=1)
                df_xvg_file.drop('solvent_restraint_energy_old', axis=1, inplace=True)
        df_xvg_file['reference'] = df_xvg_file[reference_to_column[current_stage]]
        df_xvg_file['000'] = df_xvg_file['reference']
        df_xvg_file['010'] = df_xvg_file['000'] + df_xvg_file['bonded_l']
        df_xvg_file['110'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l']
        df_xvg_file['111'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l'] + df_xvg_file['restraint_l']
        df_xvg_file['11p01'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l'] + 0.01 * df_xvg_file['restraint_l']
        df_xvg_file['11p02'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l'] + 0.02 * df_xvg_file['restraint_l']
        df_xvg_file['11p04'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l'] + 0.04 * df_xvg_file['restraint_l']
        df_xvg_file['11p08'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l'] + 0.08 * df_xvg_file['restraint_l']
        df_xvg_file['11p16'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l'] + 0.16 * df_xvg_file['restraint_l']
        df_xvg_file['11p32'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l'] + 0.32 * df_xvg_file['restraint_l']
        df_xvg_file['11p64'] = df_xvg_file['000'] + df_xvg_file['bonded_l'] + df_xvg_file['vdw_l'] + 0.64 * df_xvg_file['restraint_l']
        df_xvg_file['reference'] = df_xvg_file[reference_to_column[current_stage]]
        
        df_xvg_file['000'] = df_xvg_file['000'] - df_xvg_file['reference']
        df_xvg_file['010'] = df_xvg_file['010'] - df_xvg_file['reference']
        df_xvg_file['110'] = df_xvg_file['110'] - df_xvg_file['reference']
        df_xvg_file['111'] = df_xvg_file['111'] - df_xvg_file['reference']
        df_xvg_file['11p01'] = df_xvg_file['11p01'] - df_xvg_file['reference']
        df_xvg_file['11p02'] = df_xvg_file['11p02'] - df_xvg_file['reference']
        df_xvg_file['11p04'] = df_xvg_file['11p04'] - df_xvg_file['reference']
        df_xvg_file['11p08'] = df_xvg_file['11p08'] - df_xvg_file['reference']
        df_xvg_file['11p16'] = df_xvg_file['11p16'] - df_xvg_file['reference']
        df_xvg_file['11p32'] = df_xvg_file['11p32'] - df_xvg_file['reference']
        df_xvg_file['11p64'] = df_xvg_file['11p64'] - df_xvg_file['reference']

        water_sanity_check = pd.array(self.water_position_analysis(current_stage), dtype ="boolean")
        df_xvg_file['water_check'] = water_sanity_check
        df_xvg_file = df_xvg_file[df_xvg_file['water_check']]
        df_xvg_file['PE'] = df_xvg_file['PE'].map('{0:.2f}'.format)
        df_xvg_file['vdw_l'] = df_xvg_file['vdw_l'].map('{0:.10f}'.format)
        df_xvg_file['bonded_l'] = df_xvg_file['bonded_l'].map('{0:.7f}'.format)
        df_xvg_file['restraint_l'] = df_xvg_file['restraint_l'].map('{0:.6f}'.format)
        df_xvg_file['000'] = df_xvg_file['000'].map('{:.7f}'.format)
        df_xvg_file['010'] = df_xvg_file['010'].map('{0:.7f}'.format)
        df_xvg_file['110'] = df_xvg_file['110'].map('{0:.7f}'.format)
        df_xvg_file['11p01'] = df_xvg_file['11p01'].map('{0:.7f}'.format)
        df_xvg_file['11p02'] = df_xvg_file['11p02'].map('{0:.7f}'.format)
        df_xvg_file['11p04'] = df_xvg_file['11p04'].map('{0:.7f}'.format)
        df_xvg_file['11p08'] = df_xvg_file['11p08'].map('{0:.7f}'.format)
        df_xvg_file['11p16'] = df_xvg_file['11p16'].map('{0:.7f}'.format)
        df_xvg_file['11p32'] = df_xvg_file['11p32'].map('{0:.7f}'.format)
        df_xvg_file['11p64'] = df_xvg_file['11p64'].map('{0:.6f}'.format)
        df_xvg_file['111'] = df_xvg_file['111'].map('{0:.6f}'.format)
        df_xvg_file['pV'] = df_xvg_file['pV'].map('{0:.6f}'.format)


        if self.reading_analysis.Uscheme == "2" and float(current_stage) < 5:
            df_xvg_file = df_xvg_file[['time', 'PE', 'vdw_l', 'bonded_l', 'restraint_l', '000', '010', '110', '111', 'pV']]
        elif self.reading_analysis.Uscheme == "1" and float(current_stage) < 5:
            df_xvg_file = df_xvg_file[['time', 'PE', 'vdw_l', 'bonded_l', 'restraint_l', '000', '111', 'pV']]
        elif self.reading_analysis.Uscheme == "3" and float(current_stage) < 5:
            df_xvg_file = df_xvg_file[['time', 'PE', 'vdw_l', 'bonded_l', 'restraint_l', '000', '010', '110', '11p01', '11p02', '11p04', '11p08', '11p16', '11p32', '11p64', '111', 'pV']]
        elif float(current_stage) < 5:
            print("Unknown parameter Uscheme")
            quit()
        if self.reading_analysis.Lscheme == "1" and float(current_stage) > 4:
            df_xvg_file = df_xvg_file[['time', 'PE', 'vdw_l', 'bonded_l', 'restraint_l', '000', '010', '111', 'pV']]
        elif self.reading_analysis.Lscheme == "2" and float(current_stage) > 4:
            df_xvg_file = df_xvg_file[['time', 'PE', 'vdw_l', 'bonded_l', 'restraint_l', '000', '010', '110', '111', 'pV']]
        elif self.reading_analysis.Lscheme == "3" and float(current_stage) > 4:
            df_xvg_file = df_xvg_file[['time', 'PE', 'vdw_l', 'bonded_l', 'restraint_l', '000', '010', '110', '11p01', '11p02', '11p04', '11p08', '11p16', '11p32', '11p64', '111', 'pV']]
        elif float(current_stage) > 4:
            print("Unknown parameter Lscheme")
        
        with open (output_xvg_file, "w") as out_xvg:
            out_xvg.write("\n".join(line for line in header))
            out_xvg.write("\n")
            out_xvg.write("\n".join(line for line in legend))
            out_xvg.write("\n")
            df_xvg_file.to_csv(out_xvg, mode='a', index=False, header=False, sep=" ")
        out_xvg.close()
        return

    def engine(self):
        for stage in self.reading_analysis.stages.split(','):
            #if f"stage{stage.strip()}" in self.reading_analysis.energies['pocket_restraints']:
            #    positional_restraint_energy_for_stage = self.reading_analysis.energies['pocket_restraints'][f"stage{stage.strip()}"]
            #else:
            #    positional_restraint_energy_for_stage = None
            #    print(f"Positional restraint energy for stage {stage.strip()} doesn't exist")
            #continue
            #quit()
            #mod_file = pathlib.Path(f"{self.reading_analysis.cycle_directory}/stage{stage}/prod/dhdl_{self.reading_analysis.Uscheme}{self.reading_analysis.Lscheme}.xvg")
            self.write_xvg(current_stage=stage.strip())




