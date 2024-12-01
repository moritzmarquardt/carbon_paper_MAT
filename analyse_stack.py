import json

import matplotlib.pyplot as plt
import numpy as np
from MembraneAnalysisToolbox.DiffusionAnalysis import DiffusionAnalysis
from MembraneAnalysisToolbox.MembraneStructures import CubicMembrane, HexagonalMembrane

"""
this file is for automatic analysis of a lot of files at the same time using nohup
so ther is no plot show or plots that are for verification purposes
"""

# All hexagonal pores
hex_resname = "resname HEX and name C1"
dod_resname = "resname DOD and name C2"
hxo_resname = "resname HXO and name C1"
paths = []

# """
for nm in [2, 3, 4, 6]:
    for s in [1, 2, 3]:
        paths.append(
            [
                "/bigpool/users/ac130484/project/finished_sim/hex/poresize/"
                + str(nm)
                + "nm_NVT/simulation_"
                + str(s)
                + "/",
                180,
                hex_resname,
                dod_resname,
            ]
        )
# """

"""
for k in ["n", "y_5", "y_10", "y_15", "y_20", "y_50", "y_99"]:
    paths.append(
        [
            "/bigpool/users/ac130484/project/finished_sim/hex/hexane_dodecane/hex_18_3_2_"
            + k
            + "/",
            180,
            hex_resname,
            dod_resname,
        ]
    )
"""

"""
for k in ["n", "y_5", "y_10", "y_15", "y_20", "y_50", "y_99"]:
    paths.append(
        [
            "/bigpool/users/ac130484/project/finished_sim/hex/hexane_hexanole/hex_18_3_2_"
            + k
            + "/",
            180,
            hex_resname,
            hxo_resname,
        ]
    )
"""

"""
for s in [1, 2, 3]:
    paths.append(
        [
            "/bigpool/users/ac130484/project/cubic_box_bigger_z/hex_dod_bigger_z/sim_"
            + str(s)
            + "/",
            180 * 1.25,
            hex_resname,
            dod_resname,
        ]
    )
"""

"""
for s in [1, 2, 3]:
    paths.append(
        [
            "/bigpool/users/ac130484/project/cubic_box_hex_dod/18_2_3_n/sim_"
            + str(s)
            + "/",
            180,
            hex_resname,
            dod_resname,
        ]
    )
"""

# paths = [
#     [
#         "/bigpool/users/ac130484/project/cubic_box_hex_dod/18_2_3_n/sim_1/",
#         180,
#         hex_resname,
#         dod_resname,
#     ]
# ]
# print(paths)
print(np.shape(paths))


def analyse_resname(selector: str, short: str):
    print(f"\n{short} analysis")

    DA.calc_passagetimes(selector)
    print(f"\t{short}-passages: " + str(len(DA.passageTimes[selector])))

    DA.save_passage_times_in_ns_to_txt(selector, short + "_passagetimes_in_ns.txt")

    DA.calc_diffusion(selector)
    print(f"\t{short}-Diffusioncoefficient: " + str(DA.D[selector]).replace(".", ","))

    fig_diff = DA.plot_diffusion(selector)
    DA.save_fig_to_results(fig=fig_diff, name="diffusion_" + short)


diff_coeffs = {}

hexagonal_structure = HexagonalMembrane(
    selector="resname C",
    L=180,
)

cubic_structure = CubicMembrane(
    selector="resname C",
    cube_arrangement=(2, 2, 2),
    cube_size=90,
    pore_radius=15,
)

for o in paths:
    path = o[0]
    l = o[1]
    fir_resname = o[2]
    sec_resname = o[3]
    print("Path: " + path + "\n")

    # STEP 1: initialise the Data into the class
    DA = DiffusionAnalysis(
        topology_file=path + "topol.tpr",
        trajectory_file=path + "traj.xtc",
        results_dir=path + "analysis/",
        analysis_max_step_size_ps=200,
        verbose=False,
        membrane=hexagonal_structure,  # alternatively cubic Structure
    )

    print(DA)

    DA.find_membrane_location()
    DA.print_membrane_location()

    analyse_resname(fir_resname, fir_resname.split(" ")[1].lower())

    analyse_resname(sec_resname, sec_resname.split(" ")[1].lower())

    DA.store_results_json()
    diff_coeffs[path] = DA.D
    plt.show()
    print("\n\n\n")

print("\n\n\n RESULTS:")
print(json.dumps(diff_coeffs, indent=4))
