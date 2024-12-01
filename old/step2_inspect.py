import numpy as np
import matplotlib.pyplot as plt
import MembraneAnalysisToolbox.funcs as tfm
import MembraneAnalysisToolbox.plot as tfmp
import os
'''
STEP 2:
- inspect the received data
- Choose boundaries,
- check for timeline
- plot sample trajectories for plausibility analysis
'''

# path = '/hugepool/data/sofia_simulationen/carbon/finished_sim/hex/poresize/2nm_NVT/simulation_3/analysis/'
# print(path)

# Ask for the paths from the user
path = input("Enter the analysis folder path: ")

for res in ["hex", "dod"]:
    print("Inspecting " + res + " trajectories")
    prefix = res + "_traj_"
    x = np.load(path + prefix + "x.npy")
    y = np.load(path + prefix + "y.npy")
    z = np.load(path + prefix + "z.npy")
    timeline = np.load(path + "timeline.npy")

    print("step size in ps: " + str(timeline[2]-timeline[1]))
    print("number of timesteps: " + str(np.size(x[0, :])))
    print("sim duration (in ns if step size is in ps): " + str(np.size(x[0, :])*(timeline[2]-timeline[1])/1000))
    print("max z-value: " + str(np.max(z)) + "; min z-value: " + str(np.min(z)))
    print("number of trajs: " + str(np.size(x[:, 0])))

    fig_z_dist, ax_z_dist = plt.subplots()
    fig_z_dist.suptitle("Histogram of z", fontsize="x-large")
    ax_z_dist.hist(z.flatten(), bins=100, density=True, alpha=0.5, label=res)
    ax_z_dist.set_xlabel("z", fontsize="x-large")
    ax_z_dist.set_ylabel("Frequency", fontsize="x-large")
    ax_z_dist.legend()
    fig_z_dist.savefig(path + res + '_z_dist.png')

    fig_x_dist, ax_x_dist = plt.subplots()
    fig_x_dist.suptitle("Histogram of x", fontsize="x-large")
    ax_x_dist.hist(x.flatten(), bins=100, density=True, alpha=0.5, label=res)
    ax_x_dist.set_xlabel("x", fontsize="x-large")
    ax_x_dist.set_ylabel("Frequency", fontsize="x-large")
    ax_x_dist.legend()
    fig_x_dist.savefig(path + res + '_x_dist.png')

    fig_y_dist, ax_y_dist = plt.subplots()
    fig_y_dist.suptitle("Histogram of y", fontsize="x-large")
    ax_y_dist.hist(y.flatten(), bins=100, density=True, alpha=0.5, label=res)
    ax_y_dist.set_xlabel("y", fontsize="x-large")
    ax_y_dist.set_ylabel("Frequency", fontsize="x-large")
    ax_y_dist.legend()
    fig_y_dist.savefig(path + res + '_y_dist.png')

    rand_trajs = np.random.randint(0,str(np.size(x[:, 0])),size=(2))

    plt.figure("z_rand_trajs_" + res)
    for i in rand_trajs:
        tfmp.plot_1dtraj(z[i, :])

    plt.figure("x_rand_trajs_" + res)
    for i in rand_trajs:
        tfmp.plot_1dtraj(x[i, :])

    plt.figure("y_rand_trajs_" + res)
    for i in rand_trajs:
        tfmp.plot_1dtraj(y[i, :])
    
    plt.show()
