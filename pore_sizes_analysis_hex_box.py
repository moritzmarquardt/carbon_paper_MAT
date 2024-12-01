import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as EPSA
import numpy as np
import time
import matplotlib.pyplot as plt

# this script is hardcoded to analyse the four poresizes in the hex structure

start = time.time()

characteristics_2nm = {
    'path': "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex_NVT/2nm_NVT/",
    'y_middle': 35,
    'y_range': 10,
    'pore_size': 2
}
characteristics_3nm = {
    'path': "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex_NVT/3nm_NVT/",
    'y_middle': 30,
    'y_range': 10,
    'pore_size': 3
}
characteristics_4nm = {
    'path': "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex_NVT/4nm_NVT/",
    'y_middle': 45,
    'y_range': 10,
    'pore_size': 4
}
characteristics_6nm = {
    'path': "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex_NVT/6nm_NVT/",
    'y_middle': 50,
    'y_range': 10,
    'pore_size': 6
}

characteristics = [characteristics_2nm, characteristics_3nm, characteristics_4nm, characteristics_6nm]

for char in characteristics:
    print("Analysising Characteristics: " + str(char))
    effective_pore_sizes = []
    for i in range(3):
        Analysis = EPSA.EffectivePoreSizeAnalysis(
            topology_file = char['path'] + 'simulation_' + str(i + 1) + '/' + 'topol.tpr',
            trajectory_file = char['path'] + 'EPS_analysis/' + 'traj_simulation_' + str(i + 1) + '.xtc',
            membrane_resnames = ['C'],
            solvent_resnames = ['HEX', 'DOD'],
            y_middle = char['y_middle'],
            y_range = char['y_range'],
            verbose = True
        )
        print('Analysing simulation ' + str(i + 1) + ' in ' + char['path'])
        Analysis.analyseConstraints()
        effective_pore_size = Analysis.calculate_effective_pore_size(strategy = "intersection") / 10 #in nm
        effective_pore_sizes.append(effective_pore_size)
        print("calculated effective pore size (" + str(char['pore_size']) + "nm) is: " + str(effective_pore_size))
        print("relation between pore size and effective pore size: " + str(effective_pore_size/char['pore_size']*100) + "%")
        Analysis.plot()
        print("Time elapsed: " + str(time.time() - start))
        # plt.show()

    print("Effective pore sizes: " + str(effective_pore_sizes))
    print("Average effective pore size: " + str(sum(effective_pore_sizes) / len(effective_pore_sizes)))
    print("standard deviation: " + str(np.std(effective_pore_sizes)))
    # plt.show()

plt.show()

# Effective pore sizes: [1.8770875820363933, 1.8765139040020649, 1.8774064801837809]
# Average effective pore size: 1.8770026554074128
# standard deviation: 0.00036930787452038415

# Effective pore sizes: [2.904145761104354, 2.9042778138836893, 2.903275593340743]
# Average effective pore size: 2.9038997227762624
# standard deviation: 0.00044460667886988886

# Effective pore sizes: [3.8864302431759135, 3.8866777020792767, 3.886777204646092]
# Average effective pore size: 3.8866283833004274
# standard deviation: 0.00014587625328225327

# Effective pore sizes: [5.86730276693878, 5.866314280366277, 5.8662338323189935]
# Average effective pore size: 5.866616959874683
# standard deviation: 0.0004860497018146613