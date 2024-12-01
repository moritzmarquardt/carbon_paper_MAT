import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as EPSA
import matplotlib.pyplot as plt
import os


# Ask for the paths from the user
path = input("Enter the path to the folder containing topol.tpr and traj.xtc files: ")

# Check if topol.tpr and traj.xtc files exist in the entered path
while not os.path.exists(path + "topol.tpr") or not os.path.exists(path + "traj.xtc"):
    print("topol.tpr or traj.xtc file not found in the specified path. Please try again.")
    path = input("Enter the path to the folder containing topol.tpr and traj.xtc files: ")

resnames = input("Enter the resnames of the solvent molecules that should be analysed separated by a space (ex.: HEX DOD): ").split()
print("Entered resnames: ", resnames)

Analysis = EPSA.EffectivePoreSizeAnalysis(
            topology_file = path + 'topol.tpr',
            trajectory_file = path + 'traj.xtc',
            membrane_resnames = ['C', 'CO'],
            solvent_resnames = resnames,
            y_middle = 35,
            y_range = 10,
            verbose = True
        )

# Analysis.analyseDensity()
Analysis.analyseDensityNormalised()
#plt.show()

output_dir = os.path.join(path, "analysis_kde")
os.makedirs(output_dir, exist_ok=True)

plt.savefig(os.path.join(output_dir, 'plot.png'), bbox_inches='tight', dpi=300)

