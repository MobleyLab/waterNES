import matplotlib.pyplot as plt 
import pathlib
import os, sys
import numpy as np


def plot_dists(metadata, system_name, ligand):
    fig = plt.figure()
    fig, axs = plt.subplots(2, 4)
    fig.suptitle(f"Distance Plots | System {system_name}")

    dists = metadata

#    for i in range(1, 3):
#        for j in range(1, 4):
    if 'stage1' in dists:
        axs[0, 0].plot(dists['stage1']['trapped water'])
        axs[0, 0].plot(dists['stage1']['solvent water 1st'])
        axs[0, 0].plot(dists['stage1']['solvent water 2nd'])
        axs[0, 0].set_title('Stage 1')

    if 'stage2' in dists:
        axs[0, 1].plot(dists['stage2']['trapped water'])
        axs[0, 1].plot(dists['stage2']['solvent water 1st'])
        axs[0, 1].plot(dists['stage2']['solvent water 2nd'])
        axs[0, 1].set_title('Stage 2')

    if 'stage3' in dists:
        axs[0, 2].plot(dists['stage3']['trapped water'])
        axs[0, 2].plot(dists['stage3']['solvent water 1st'])
        axs[0, 2].plot(dists['stage3']['solvent water 2nd'])
        axs[0, 2].set_title('Stage 3')

    if 'stage4' in dists:
        axs[0, 3].plot(dists['stage4']['trapped water'])
        axs[0, 3].plot(dists['stage4']['solvent water 1st'])
        axs[0, 3].plot(dists['stage4']['solvent water 2nd'])
        axs[0, 3].set_title('Stage 4')

    if 'stage5' in dists:
        axs[1, 0].plot(dists['stage5']['trapped water'])
        axs[1, 0].plot(dists['stage5']['solvent water 1st'])
        axs[1, 0].plot(dists['stage5']['solvent water 2nd'])
        axs[1, 0].set_title('Stage 5')

    if 'stage6' in dists:
        axs[1, 1].plot(dists['stage6']['trapped water'])
        axs[1, 1].plot(dists['stage6']['solvent water 1st'])
        axs[1, 1].plot(dists['stage6']['solvent water 2nd'])
        axs[1, 1].set_title('Stage 6')

    if 'stage7' in dists:
        axs[1, 2].plot(dists['stage7']['trapped water'])
        axs[1, 2].plot(dists['stage7']['solvent water 1st'])
        axs[1, 2].plot(dists['stage7']['solvent water 2nd'])
        axs[1, 2].set_title('Stage 7')

    
    for ax in axs.flat:
        ax.set(xlabel='time (ps)', ylabel='distance (nm)', xticks=np.arange(0, 700, 200), yticks=np.arange(0, 10, 1))
#        ax.tick_params(axis='x', nticks=5)
        ax.label_outer()
    pathlib.Path(os.path.join("/dfs9/dmobley-lab/swapnilw/openeye_colab/results", "_".join(system_name.split('_')[:2]))).mkdir(parents=True, exist_ok=True)
    out_file_path = os.path.join("/dfs9/dmobley-lab/swapnilw/openeye_colab/results", "_".join(system_name.split('_')[:2]), f"distances_{system_name.split('_')[0]}_{system_name.split('_')[3]}_1-7.png")
    if ligand == "A":
        pass
    else:
        out_file_path = os.path.join("/dfs9/dmobley-lab/swapnilw/openeye_colab/results", "_".join(system_name.split('_')[:2]), f"distances_{system_name.split('_')[1]}_{system_name.split('_')[3]}_1-7.png")
    plt.savefig(out_file_path)
    
    fig.clear()

    if 'stage6.1' in dists:
        fig, axs = plt.subplots(2, 4)
        fig.suptitle(f"Distance Plots | System {system_name}")
        axs[0, 0].plot(dists['stage6.1']['trapped water'])
        axs[0, 0].plot(dists['stage6.1']['solvent water 1st'])
        axs[0, 0].plot(dists['stage6.1']['solvent water 2nd'])
        axs[0, 0].set_title('Stage 6.1')
    
        axs[0, 1].plot(dists['stage6.2']['trapped water'])
        axs[0, 1].plot(dists['stage6.2']['solvent water 1st'])
        axs[0, 1].plot(dists['stage6.2']['solvent water 2nd'])
        axs[0, 1].set_title('Stage 6.2')
    
        axs[0, 2].plot(dists['stage6.3']['trapped water'])
        axs[0, 2].plot(dists['stage6.3']['solvent water 1st'])
        axs[0, 2].plot(dists['stage6.3']['solvent water 2nd'])
        axs[0, 2].set_title('Stage 6.3')
    
        axs[0, 3].plot(dists['stage6.4']['trapped water'])
        axs[0, 3].plot(dists['stage6.4']['solvent water 1st'])
        axs[0, 3].plot(dists['stage6.4']['solvent water 2nd'])
        axs[0, 3].set_title('Stage 6.4')
    
        axs[1, 0].plot(dists['stage6.5']['trapped water'])
        axs[1, 0].plot(dists['stage6.5']['solvent water 1st'])
        axs[1, 0].plot(dists['stage6.5']['solvent water 2nd'])
        axs[1, 0].set_title('Stage 6.5')
    
        axs[1, 1].plot(dists['stage6.6']['trapped water'])
        axs[1, 1].plot(dists['stage6.6']['solvent water 1st'])
        axs[1, 1].plot(dists['stage6.6']['solvent water 2nd'])
        axs[1, 1].set_title('Stage 6.6')
    
        axs[1, 2].plot(dists['stage6.7']['trapped water'])
        axs[1, 2].plot(dists['stage6.7']['solvent water 1st'])
        axs[1, 2].plot(dists['stage6.7']['solvent water 2nd'])
        axs[1, 2].set_title('Stage 6.7')

        for ax in axs.flat:
            ax.set(xlabel='time (ps)', ylabel='distance (nm)', xticks=np.arange(0, 700, 200), yticks=np.arange(0, 10, 1))
#            ax.tick_params(axis='x', nticks=5)
            ax.label_outer()
            
        out_file_path = os.path.join("/dfs9/dmobley-lab/swapnilw/openeye_colab/results", system_name.split('_')[0], f"distances_{system_name.split('_')[1]}_6.1-6.7.png")
        plt.savefig(out_file_path)
        fig.clear()



