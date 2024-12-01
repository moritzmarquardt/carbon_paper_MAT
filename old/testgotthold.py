import matplotlib.pyplot as plt
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.optimize import least_squares


# path = 
# prefix = 
# ffs = np.load(path + prefix + "ffs_transition.npy")
# ffe = np.load(path + prefix + "ffe_transition.npy")

series = []

with open('/bigpool/users/st166545/TransitionAnalysisProject/tt-ns_hex-dod_L181_310K_N283.txt', 'r') as file:
    for line in file:
        series.append(float(line.strip()))

series = list(series) # ? neccesary?

print(len(series))

L = 181 #length of the path that the transitions took
N = 0 # not important for diffusion
T = 310 # in kelvin 

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

def hom(x,D,i,L):
    t=(L)**2/(i**2*np.pi**2*D)
    return((-1)**(i-1)*i**2*np.exp(-x/t))

def fitfunc_hom(x,D,L):
    i=151
    result=0
    for j in range(1,i):
        result=result+hom(x,D,j,L)
    return(2*np.pi**2*D/(L)**2*result)

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