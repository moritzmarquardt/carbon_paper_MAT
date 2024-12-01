import matplotlib.pyplot as plt
import numpy as np
from MembraneAnalysisToolbox.DiffusionAnalysis import DiffusionAnalysis

"""
this file is for testing and development purposes only of the bootstrapping function in the TPAnalysis class
"""
# paths = [str(nm) + "nm_NVT/simulation_" + str(s) + "/" for nm in [2,3,4,6] for s in [1,2,3]]
# paths = ["3nm_NVT/simulation_2/"]
paths = [
    str(nm) + "nm_NVT/simulation_" + str(s) + "/"
    for nm in [2, 3, 4, 6]
    for s in [1, 2, 3]
]

for path in paths:
    path = "/bigpool/users/ac130484/project/finished_sim/hex/poresize/" + path
    print("Path: " + path)

    # STEP 1: initialise the Data into the class
    DA = DiffusionAnalysis(
        topology_file=path + "topol.tpr",
        trajectory_file=path + "traj.xtc",
        analysed_max_step_size_ps=200,
        results_dir=path + "analysis/",
        verbose=True,
        L=180,
        # z_lower=233.23501586914062,
        # D=16.709077586864158,
    )
    print(DA)

    print("Analysis of the membrane z-dimension")
    DA.find_membrane_location_hexstructure(mem_selector="resname C")
    print("\tz_lower: " + str(DA.z_lower))

    print("\nHEX analysis")

    ffs, ffe = DA.calc_passagetimes(["resname HEX and name C1"])
    print("\tpassages: " + str(len(ffs)))
    D_hex = DA.calc_diffusion(list(ffe - ffs))
    print("\tDiffusioncoefficient: " + str(D_hex).replace(".", ","))

    # BOOTSTRAPPING Kri
    bootstrap_diffs = DA.bootstrap_diffusion(
        "resname HEX and name C1", n_bootstraps=2, plot=False
    )

    # Bootstrapping internet
    # bootstrap_diffs = Analysis2nm_1.bootstrapping_diffusion(
    #     selector="resname HEX and name C1",
    #     bootstrap_sample_length_ns=10,
    #     n_bootstraps=1000,
    #     z_lower=z_lower,
    #     L=L,
    #     plot=False,
    # )

    print("\nBootstraped Diffusion Coefficients: " + str(bootstrap_diffs))

    plt.hist(bootstrap_diffs)
    plt.axvline(D_hex, color="r", linestyle="dashed", linewidth=1)

    print("mean: " + str(np.mean(bootstrap_diffs)))
    print("std: " + str(np.std(bootstrap_diffs)))
    print(
        "standard error of the mean: "
        + str(np.std(bootstrap_diffs) / np.sqrt(len(bootstrap_diffs)))
    )

    plt.show()
