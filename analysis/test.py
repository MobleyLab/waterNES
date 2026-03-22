import os, sys
import pickle
sys.path.append('/dfs4/dmobley-lab/swapnilw/waterNES/')
from water_nes.analysis.free_energy_estimate import FreeEnergyEstimate
from water_nes.analysis.nes_free_energy import calculate_nes_free_energy


with open ('analysistest_dist.pickle', 'rb') as f:
        analysis = pickle.load(f)
f.close()
breakpoint()
print(analysis)
