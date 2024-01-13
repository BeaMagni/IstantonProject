import numpy as np
from numba import njit
from scipy.integrate import simps
import General_functions as fn

etha = float(input("Insert the value for etha, the shift of the potential: "))
w0 = 4*etha
a = float(input("Insert the value for a, the lattice spacing: "))
delta_x = 0.5
n_equil = int(input("Insert the number of Monte Carlo equilibration sweeps: "))
n_sweeps = int(input("Insert the number of Monte Carlo sweeps (must be bigger than the previous value): "))
n_switching = int(input("Insert the number of adiabatic switching sweeps: "))
n_beta = int(input("Insert the number of steps for beta: "))
beta_max = float(input("Insert the maximum value for beta: "))
beta = np.linspace(1.0,beta_max,n_beta)
temperature = 1.0/beta
F = np.zeros(n_beta)
F_err = np.zeros(n_beta)
start = bool(input("Insert the desired start (0-cold, 1-hot): "))

for b in range(n_beta):
    N = int(beta[b]/a)
    F[b], F_err[b] = fn.montecarlo_switching(N,n_equil,n_sweeps,n_switching,etha,start,a,delta_x)
    F[b] /= beta[b]
    F[b] += fn.free_energy_zero(beta[b],w0)
    F_err[b] /= beta[b]

np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/adiabatic/temperature.txt",temperature)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/adiabatic/free_energy.txt",F)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/adiabatic/free_energy_err.txt",F_err)
