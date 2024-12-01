import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as mat

'''
this file is for testing and development purposes only
'''

path = "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/hex/3nm_1musonly/simulation_1/" # diff pore sizes directory, ymiddle 40
# path = "/bigpool/data/projects/Carbon_pores_Sofia/correctdensity/cubic/3nm/simulation_1/" 

Analysis3nm = mat.EffectivePoreSizeAnalysis( 
    topology_file = path + 'topol.tpr', 
    trajectory_file = path + 'traj.xtc',
    membrane_resnames = ['C'],
    solvent_resnames = ['HEX', 'DOD'],
    y_middle = 45, 
    y_range = 10, 
    verbose = True
)
Analysis3nm.analyseConstraints()
print("calculated effective pore size (3nm) is: " + str(Analysis3nm.calculate_effective_pore_size(strategy = "intersection") / 10) + "nm")
Analysis3nm.plot()