# inspct file that has routines and workflow to inspect trajectories
# shoukd be possible too read out boundaries,
# verify data
import numpy as np
import matplotlib.pyplot as plt
import MembraneAnalysisToolbox.funcs as tfm
import MembraneAnalysisToolbox.plot as tfmp

# path = "/hugepool/data/sofia_simulationen/carbon/finished_sim/hex/poresize/2nm_NVT/simulation_3/analysis/"
# print(path)

# Ask for the paths from the user
path = input("Enter the analysis folder path: ")

for res in ["hex", "dod"]:
    prefix = res + "_traj_"
    indexes = np.load(path + prefix + "indizes_transition.npy")
    x_passages = np.load(path + prefix + "x.npy")[indexes]
    y_passages = np.load(path + prefix + "y.npy")[indexes]
    z_passages = np.load(path + prefix + "z.npy")[indexes]
    ffs = np.load(path + prefix + "ffs_transition.npy")
    ffe = np.load(path + prefix + "ffe_transition.npy")
    # distances = np.load(path + prefix + "distances_hor.npy") 
    timeline = np.load(path + "timeline.npy")

    plt.figure("Verteilung der Durchgangszeiten")
    tfmp.plot_dist(ffe-ffs,number_of_bins=10, max_range=np.max(ffe-ffs))
    plt.xlabel("Durchgangszeiten")
    plt.ylabel("relative Häufigkeit")

    #plot 3d trajectories
    kernel_size = 1000
    kernel = np.ones(kernel_size) / kernel_size
    fig = plt.figure("3d trajektorien")
    ax = fig.add_subplot(projection='3d')
    fig2, ax2 = plt.subplots()
    for sel in np.random.randint(0,x_passages.shape[0],size=(3)):
        print(sel)
        ax2.plot(np.arange(z_passages[sel,ffs[sel]-1:ffe[sel]+2].size),z_passages[sel,ffs[sel]-1:ffe[sel]+2], label=str(sel))
        x_passages_sel = x_passages[sel]
        y_passages_sel = y_passages[sel]
        z_passages_sel = z_passages[sel]
        slicer = slice(ffs[sel]-1, ffe[sel]+2, 1)
        tfmp.plot_3dtraj(ax, x_passages_sel[slicer],y_passages_sel[slicer],z_passages_sel[slicer]) #++1 so that the last timestep is also included
        tfmp.plot_3dpoints(ax, x_passages_sel[ffs[sel]-1],y_passages_sel[ffs[sel]-1],z_passages_sel[ffs[sel]-1]) #starting point
        tfmp.plot_3dpoints(ax, x_passages_sel[ffe[sel]+2],y_passages_sel[ffe[sel]+2],z_passages_sel[ffe[sel]+2]) #end point
    # tfmp.plot_3dbounds(ax, zbounds)


    # plot starting points
    fig = plt.figure("plotten aller Startpunkte")
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel("x in nm", fontsize="x-large")
    ax.set_ylabel("y in nm", fontsize="x-large")
    ax.set_zlabel("z in nm", fontsize="x-large")
    ax.set_title("Membrane-entry points of the passage-trajectories (" + res.upper() + ")", fontsize="x-large")
    ax.scatter(x_passages[np.arange(np.size(x_passages,0)),ffs+1]/10,y_passages[np.arange(np.size(x_passages,0)),ffs+1]/10,z_passages[np.arange(np.size(x_passages,0)),ffs+1]/10) #ugly way of getting the point. maybe there is a better way
    # tfmp.plot_3dpoints(ax,x_passages[np.arange(np.size(x_passages,0)),ffs+200],y_passages[np.arange(np.size(x_passages,0)),ffs+200],z_passages[np.arange(np.size(x_passages,0)),ffs+200]) #ugly way of getting the point. maybe there is a better way
    # tfmp.plot_3dpoints(ax,x_passages[np.arange(np.size(x_passages,0)),ffs+3000],y_passages[np.arange(np.size(x_passages,0)),ffs+3000],z_passages[np.arange(np.size(x_passages,0)),ffs+3000]) #ugly way of getting the point. maybe there is a better way

    plt.show()




#Verteilung der horizontal zurückgelegten Strecke
# plt.figure("Verteilung der quer zurückgelegten Strecke")
# tfmp.plot_dist(distances,number_of_bins=30,max_range=np.max(distances))
# plt.xlabel("horizontale strecke dodecane")
# plt.ylabel("relative Häufigkeit")

# direct = tfm.path_cat(x_passages,y_passages,ffs,ffe)
# print("direkt Durchgänge: " + str(direct))


#plotten einer 3d trajectorie
# sel = 1



