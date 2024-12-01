import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as EPSA
import matplotlib.pyplot as plt


for res in ['HEX', 'HXO']:
    Analysis = EPSA.EffectivePoreSizeAnalysis(
                topology_file = '/bigpool/users/ac130484/project/hex_box_hex_hexanol/hex_18_2_3_y_99/topol.tpr',
                trajectory_file = '/bigpool/users/ac130484/project/hex_box_hex_hexanol/hex_18_2_3_y_99/traj.xtc',
                membrane_resnames = ['C'],
                solvent_resnames = [res],
                y_middle = 30,
                y_range = 10,
                verbose = True
            )

    # Analysis.analyseConstraints()
    # plt.show()
    Analysis.analyseDensityNormalised(slice_height = 10, factors = ['scott', 'scott'], skip = 5)
plt.show()