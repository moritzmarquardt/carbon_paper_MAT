import numpy as np
import MembraneAnalysisToolbox.funcs as tfm
import matplotlib.pyplot as plt

# Diese Datei ist dafür da Stichprobenmäßig zu überprüfen, ob die funktion dur dist improved funktioniert

# Ask for the paths from the user
path = "/hugepool/data/sofia_simulationen/carbon/finished_sim/hex/poresize/2nm_NVT/simulation_3/analysis/"

lower_zbound = 240
upper_zbound = 400


timeline = np.load(path + "timeline.npy")
zbounds = [lower_zbound, upper_zbound]  # obtained with the inspect data script

res = "hex"
prefix = res + "_traj_"

z = np.load(path + prefix + "z.npy")

# GET PASSAGES AND TRANSITION DURATION

random_transition = np.random.randint(0, z.shape[0])
# random_transition = 503
print("Random transition: ", random_transition)
plt.plot(z[random_transition].T)
plt.plot(z[random_transition+1].T)
plt.show()
ffs, ffe, indizes = tfm.dur_dist_improved(z[random_transition:random_transition+2], zbounds)
print(ffs, ffe, indizes)
# plt.plot(z[random_transition:random_transition+1, :].T)
# plt.scatter(ffs, z[random_transition, ffs], c='r')
# plt.scatter(ffe, z[random_transition, ffe], c='g')
