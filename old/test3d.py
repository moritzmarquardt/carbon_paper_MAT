import numpy as np
import matplotlib.pyplot as plt
import MembraneAnalysisToolbox.funcs as tfm
import MembraneAnalysisToolbox.plot as tfmp
import time

start = time.time()

wsl_path = '/home/moritz/Projektarbeit/3dtrajs_zusatz/'  #paths for wsl vs code enviroment (when you open the foolder/workspace in wsl mode)
windows_path = '//wsl.localhost/Ubuntu/home/moritz/Projektarbeit/3dtrajs_zusatz/' #paths for default windows vs coode enviroment
root_path = wsl_path
postfix = '_nojump_all'   #postfix als dateiendung: _all sind alle trajs, _nojump sind die trajs ohne pbc jumps, '' ist random trajs
print(postfix)

#path building
x_dod = np.load(root_path + 'testdata/x_dod_c2' + postfix + '.npy')
y_dod = np.load(root_path + 'testdata/y_dod_c2' + postfix + '.npy')
z_dod = np.load(root_path + 'testdata/z_dod_c2' + postfix + '.npy')
x_hex = np.load(root_path + 'testdata/x_hex_c1' + postfix + '.npy')
y_hex = np.load(root_path + 'testdata/y_hex_c1' + postfix + '.npy')
z_hex = np.load(root_path + 'testdata/z_hex_c1' + postfix + '.npy')
timeline = np.load(root_path + 'testdata/timeline.npy')

number_of_times = timeline.size
zbounds = [20.5,38.5] #von sofia
xbounds =  [1.5,3,6,7.5,10.5,12,15,16.5] #aus eigenen berechnungen und sofia rücksprache


#dod
ffs, ffe, indizes = tfm.dur_dist_improved(z_dod,zbounds,p_middle=1)
print("passages: " + str(ffs.size))
x_dod_passages = x_dod[indizes]
y_dod_passages = y_dod[indizes]
z_dod_passages = z_dod[indizes]

#Z-traj plotten
'''sel = 10
plt.figure("Z-trajectorie des " + str(sel) + ". Durchgangs")
tfmp.plot_1dtraj(z_dod_passages[sel])
tfmp.plot_hor_bounds(zbounds)
tfmp.plot_point(ffs[sel],z_dod_passages[sel,ffs[sel]])
tfmp.plot_point(ffe[sel],z_dod_passages[sel,ffe[sel]])
plt.xlabel("Zeitschritte")
plt.ylabel("Z-Wert")'''

#plotten einer 3d trajectorie
'''sel = 10
fig = plt.figure("3d trajektorie des " + str(sel) + ". Durchgangs")
ax = fig.add_subplot(projection='3d')
tfmp.plot_3dtraj(ax, x_dod_passages[sel][ffs[sel]:ffe[sel]+1],y_dod_passages[sel][ffs[sel]:ffe[sel]+1],z_dod_passages[sel][ffs[sel]:ffe[sel]+1]) #++1 so that the last timestep is also included
tfmp.plot_3dpoints(ax, x_dod_passages[sel,ffs[sel]],y_dod_passages[sel,ffs[sel]],z_dod_passages[sel,ffs[sel]]) #starting point
tfmp.plot_3dpoints(ax, x_dod_passages[sel,ffe[sel]],y_dod_passages[sel,ffe[sel]],z_dod_passages[sel,ffe[sel]]) #end point
tfmp.plot_3dbounds(ax, xbounds)'''


#Durchgangszeten verteilung plotten
'''plt.figure("Verteilung der Durchgangszeiten")
tfmp.plot_dist(ffe-ffs,number_of_bins=30, max_range=1000)
plt.xlabel("Durchgangszeiten")
plt.ylabel("relative Häufigkeit")'''


distances = tfm.calc_hor_dist(x_dod_passages,y_dod_passages,ffs,ffe)
print("Gesamte Distanz: " + str(np.sum(distances)))

#Verteilung der horizontal zurückgelegten Strecke
plt.figure("dodecane Verteilung der quer zurückgelegten Strecke")
tfmp.plot_dist(distances,number_of_bins=30,max_range=800)
plt.xlabel("horizontale strecke dodecane")
plt.ylabel("relative Häufigkeit")

# plot starting points
'''fig = plt.figure("plotten aller Startpunkte")
ax = fig.add_subplot(projection='3d')
tfmp.plot_3dpoints(ax,x_dod_passages[np.arange(np.size(x_dod_passages,0)),ffs+1],y_dod_passages[np.arange(np.size(x_dod_passages,0)),ffs+1],z_dod_passages[np.arange(np.size(x_dod_passages,0)),ffs+1]) #ugly way of getting the point. maybe there is a better way
'''

direct = tfm.path_cat(x_dod_passages,y_dod_passages,ffs,ffe)
print("direkt Durchgänge: " + str(direct))
print("elapsed time: " + str(time.time() - start))
plt.show()