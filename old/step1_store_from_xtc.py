import numpy as np
import MDAnalysis as mda
import os
import json

'''
STEP 1:
- store the numpy files of x,y,z trajectories for both solvents from a traj.xtc file
- it is assumed that the traj.xtc file is in the same folder as the topol.tpr file
- The files (that end with xvg - so trajectory files) are stored into the same folder with the same name but with the ending .npy
'''


# path = '/bigpool/users/ac130484/project/finished_sim/hex/poresize/2nm_NVT/simulation_3/'
# save_path = '/hugepool/data/sofia_simulationen/carbon/finished_sim/hex/poresize/2nm_NVT/simulation_3/analysis/'

# Ask for the paths from the user
path = input("Enter the path to the folder containing topol.tpr and traj.xtc files: ")

# Check if topol.tpr and traj.xtc files exist in the entered path
while not os.path.exists(path + "topol.tpr") or not os.path.exists(path + "traj.xtc"):
    print("topol.tpr or traj.xtc file not found in the specified path. Please try again.")
    path = input("Enter the path to the folder containing topol.tpr and traj.xtc files: ")

save_path = input("Enter the path to the folder where the analysis and .npy files will be saved: ")

# Check if topol.tpr and traj.xtc files exist in the path
if not os.path.exists(path + "topol.tpr") or not os.path.exists(path + "traj.xtc"):
    raise Exception("topol.tpr or traj.xtc file not found in the specified path -> storing is not possible")

print("store traj.xtc file as .npy files")
# the topol file contains the information about the system
# the traj file contains the trajectory of the system
u = mda.Universe(path + "topol.tpr", path + 'traj.xtc')
print("number of frames: " + str(u.trajectory.n_frames))
print("step size: " + str(u.trajectory.dt) + " ps")
print("simulation duration: " + str((u.trajectory.n_frames-1)*u.trajectory.dt/1000) + " ns")

# select c2 atoms of the hex and dod molecules because otherwise 
# the passages of each oof the molecuele atoms are being counted
hex = u.select_atoms("resname HEX and name C1") # hex has 2 beats per molecue in the coarse grained model, sofia always used c1 for hex
dod = u.select_atoms("resname DOD and name C2") # dod has 3 beats per molecue in the coarse grained model so select the middle one
print("number of hex atoms: " + str(hex.n_atoms))
print("number of dod atoms: " + str(dod.n_atoms))

#we want at least 5 steps per ns because transitions can last only 1ns and it needs steps to identify it as a passage
nth = int(np.floor(1000/u.trajectory.dt/5))
print("every " + str(nth) + "th frame is stored")
timesteps = int(np.ceil(u.trajectory.n_frames/nth))
hex_traj = np.zeros((hex.n_atoms, timesteps, 3))
dod_traj = np.zeros((dod.n_atoms, timesteps, 3))
timeline = np.zeros(timesteps)
for i,ts in enumerate(u.trajectory[::nth]):
    hex_traj[:,i,:] = hex.positions
    dod_traj[:,i,:] = dod.positions
    timeline[i] = ts.time

for i,j in enumerate(["x", "y", "z"]):
    print("save " + j + " trajectories")
    np.save(
        save_path + "hex_traj_" + j + ".npy",
        hex_traj[:,:,i],
    )
    print(hex_traj[:,:,i].shape)
    np.save(
        save_path + "dod_traj_" + j + ".npy",
        dod_traj[:,:,i],
    )
    print(dod_traj[:,:,i].shape)

results = {
    "topology_path": path + "topol.tpr",
    "trajectory_path": path + "traj.xtc",
    "output_path": save_path,
    "number_of_frames": u.trajectory.n_frames,
    "step_size": u.trajectory.dt,
    "simulation_duration": (u.trajectory.n_frames-1)*u.trajectory.dt/1000,
    "number_of_hex_atoms": hex.n_atoms,
    "number_of_dod_atoms": dod.n_atoms,
    "nth": nth,
}
    
np.save(save_path + "analysis_results.json", results)
np.save(save_path + "timeline.npy", timeline)
print("saving is done")