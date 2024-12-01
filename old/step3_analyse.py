import numpy as np
import matplotlib.pyplot as plt
import MembraneAnalysisToolbox.funcs as tfm
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.optimize import least_squares


# CDF #####################################################

def hom_cdf(x,D,i,L):
    t=(L)**2/(i**2*np.pi**2*D) #L^2/(i^2*pi^2*D)
    return((-1)**(i-1)*np.exp(-x/t)) #summand in Gl. 10 vanHijkoop

def fitfunc_hom_cdf(x,D,L):
    i=50 #Summe geht bis 50 (approx statt undendlich)
    result=0
    for j in range(1,i):
        result=result+hom_cdf(x,D,j,L)
    return(1-2*result) #gleichung 10 in vanHijkoop paper

def fitfunc_hom_cdf_lsq(L):
    def f(D,x,y):
        i=50
        result=0
        for j in range(1,i):
            result=result+hom_cdf(x,D,j,L)
        return(1-2*result-y)
    return f

def fitting_hom_cdf_lsq(x_data,y_data, L):
    res_robust = least_squares(fitfunc_hom_cdf_lsq(L), x0=20, loss ='soft_l1', f_scale=0.3, args=(x_data, y_data))
    return res_robust.x


# PDF #####################################################

def hom(x,D,i,L):
    t=(L)**2/(i**2*np.pi**2*D)
    return((-1)**(i-1)*i**2*np.exp(-x/t))

def fitfunc_hom(x,D,L):
    i=151
    result=0
    for j in range(1,i):
        result=result+hom(x,D,j,L)
    return(2*np.pi**2*D/(L)**2*result)


# path = '/hugepool/data/sofia_simulationen/carbon/finished_sim/hex/poresize/2nm_NVT/simulation_3/analysis/'

# Ask for the paths from the user
path = input("Enter the analysis folder path: ")

lower_zbound = float(input("Enter the lower z-boundary: "))
upper_zbound = float(input("Enter the upper z-boundary: "))


timeline = np.load(path + "timeline.npy")
zbounds = [lower_zbound, upper_zbound]  # obtained with the inspect data script

for res in ["hex", "dod"]:
    prefix = res + "_traj_"

    x = np.load(path + prefix + "x.npy")
    y = np.load(path + prefix + "y.npy")
    z = np.load(path + prefix + "z.npy")

    # GET PASSAGES AND TRANSITION DURATION

    ffs, ffe, indizes = tfm.dur_dist_improved(z, zbounds)
    np.save(path + prefix + "indizes_transition.npy",indizes)
    np.save(path + prefix + "ffs_transition.npy",ffs)
    np.save(path + prefix + "ffe_transition.npy",ffe)
    # print(ffs)
    # print(ffe)
    # print(indizes)
    print(res + "-passages: " + str(ffs.size)) 

    series = list(ffe - ffs)

    L = upper_zbound - lower_zbound #length of the path that the transitions took
    T = 296 # in kelvin 

    # CDF #####################################################

    ecdf = ECDF(series)
    idx = (np.abs(ecdf.y - 0.5)).argmin()
    centertime = ecdf.x[idx]

    """ FIT DATA """
    params_hom_cdf = fitting_hom_cdf_lsq(ecdf.x[1:],ecdf.y[1:], L)

    """ PLOT DATA """
    x_lim = centertime*4
    x_cdf = np.linspace(0,x_lim*2,200)
    y_hom_cdf = fitfunc_hom_cdf(x_cdf, params_hom_cdf[0], L)

    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle('PDF and CDF fit')
    ax2.scatter(ecdf.x, ecdf.y, color=[0,0.5,0.5])
    ax2.plot(x_cdf[1:],y_hom_cdf[1:],label='hom', color='red', ls = 'dashed')
    ax2.legend(loc='center right')
    ax2.set_xlim(0,x_lim)

    D_hom_cdf=params_hom_cdf[0]
    print("D_hom_cdf: " + str(D_hom_cdf))

    # PDF #####################################################

    """ PREPARE DATA """
    bins = int(10*np.max(series)/centertime)
    histo, edges = np.histogram(series, bins, density=True);
    center=edges-(edges[2]-edges[1]);
    center=np.delete(center,0)
    edges=np.delete(edges,0)

    params_hom = params_hom_cdf

    """ PLOT DATA """
    x = x_cdf
    y2 = fitfunc_hom(x,params_hom[0], L)

    ax1.hist(series,bins=len(center), density=True, color=[0,0.5,0.5])
    ax1.plot(x[1:],y2[1:],label='hom', color='red', ls='dashed')
    ax1.set_xlim(0,x_lim)
    ax1.set_ylim(0,1.2*np.max(histo[0:]))
    ax1.legend(loc='center right')

    D_hom=params_hom[0]
    print("D_hom: " + str(D_hom))

    plt.show()


    # distances = tfm.calc_hor_dist(x_passages,y_passages,ffs,ffe)
    # print("Gesamte Distanz: " + str(np.sum(distances)))
    # np.save(path + prefix + "distances_hor.npy", distances)
