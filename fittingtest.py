import MembraneAnalysisToolbox.EffectivePoreSizeAnalysis as mat
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import scipy.misc as spm

def integrand(x,R,H):
    return H - np.sqrt(R**2 - x**2)

def func(x, R, H):
    # return 2 * spi.quad(integrand, -R, x, args=(R, H))[0] #numerical integration
    def bfunc(x,R,H):
        return -R**2*np.arcsin(x/R)-x*np.sqrt(R**2-x**2) + 2*H*x
    return  (bfunc(x,R,H) - bfunc(-R,R,H)) # analytical integration


def cum(x, R, H, a):
    """
    Calculate the cumulative value based on the given parameters.

    Parameters:
    x (float): The input value.
    R (float): The radius of the pore
    H (float): 
    a (float): The offset value.

    Returns:
    float: The calculated cumulative value.
    """
    if x < -R:
        return H * (x + a)
    elif x > R:
        return H * (a - R) + 2 * R * 2 * H - np.pi * R**2 + H * (x - R)
    elif x <= R and x >= -R:
        return func(x, R, H) + H * (a - R)
    

R = 1
H = 1
a = 2
x = np.linspace(-a, a, 100)
plt.plot(x, [cum(val, R, H, a) for val in x])
# derivative = [spm.derivative(cum, val, dx=1e-6, args=(R, H, a)) for val in x]
# plt.plot(x, derivative)
# plt.plot(x, [integrand(val, R, H) for val in x])
plt.show()
