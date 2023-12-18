import numpy as np
from numba import njit
from scipy.integrate import simps

@njit
def initialize_lattice(n,etha,start):
    if start is True:
        x = np.zeros(n)
        for j in range(n):
            x[j] = -etha
        return x
    elif start is False:
        x = np.random.uniform(-etha,etha,size = n)
        x[0] = x[-2]
        x[-1] = x[1]
        return x

@njit
def potential_harmonic(x,w0):
    potential = (1/4.0)*np.power(x*w0,2)
    return potential

@njit
def potential_anharmonic(x, etha):
    potential = np.power(np.power(x,2)-np.power(etha,2),2)
    return potential

@njit
def stat_av_var(observable,observable2,n_data):
    observable_av = observable/n_data
    observable_err = np.sqrt(observable2/np.power(n_data,2)-np.power(observable_av,2)/n_data)
    return observable_av,observable_err
