import numpy as np
from scipy.integrate import simps #we use this integration function since it gives better results than others
import General_functions as fn

#monte carlo procedure in the adiabatic switching case

def montecarlo_switching(n, n_equil, n_sweeps, n_switching, etha, start, a, delta_x):
    w0 = 4*etha
    d_alpha = 1.0/n_switching
    delta_s_alpha = np.zeros(2*n_switching+1)
    delta_s_alpha2 = np.zeros(2*n_switching+1)
    x = fn.lattice_initialization(n,etha,start)
    
    for i_switching in range(2*n_switching+1):
        if i_switching <= n_switching:
            alpha = i_switching*d_alpha
        else:
            alpha = 2.0 - i_switching*d_alpha
        
        for _ in range(n_equil): #equilibration runs for the metropolis algorithm
            fn.metropolis_switching(x,etha,a,delta_x,alpha)
        
        for _ in range(n_sweeps-n_equil):
            delta_s_alpha_temp = 0.0
            fn.metropolis_switching(x,etha,a,delta_x,alpha)
            potential_diff = fn.potential_anharmonic(x,etha)
            potential_diff -= fn.potential_harmonic(x,w0)
            delta_s_alpha_temp = np.sum(potential_diff)*a
            delta_s_alpha[i_switching] += delta_s_alpha_temp
            delta_s_alpha2[i_switching] += np.power(delta_s_alpha_temp,2)
    
    delta_s_alpha_av, delta_s_alpha_err = fn.average_std(delta_s_alpha,delta_s_alpha2,n_sweeps-n_equil)
    integral = simps(delta_s_alpha_av[0:n_switching+1],dx=d_alpha)
    integral_inv = simps(delta_s_alpha_av[n_switching:2*n_switching+1],dx=d_alpha) #we take into account the hyseresis effects to reduce the presence of systematics errors
    integral_err = simps(delta_s_alpha_err[0:n_switching+1],dx=d_alpha)
    integral_inv_err = simps(delta_s_alpha_err[n_switching:2*n_switching+1],dx=d_alpha)
    propagation_error = np.sqrt(integral_err+integral_inv_err)/2.0
    hysteresis_error = np.abs(integral-integral_inv) 
    result = -(integral+integral_inv)/2.0
    result_err = np.sqrt(np.power(propagation_error,2)+np.power(hysteresis_error,2))
    return result, result_err
    
#explicit computation of the free energy using tha adiabatic switching procedure. The parameters are left as an input from keyboard.

def main():
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
    
    output_path = './instanton_project/adiabatic'
    fn.path_creation(output_path)

    for b in range(n_beta):
        N = int(beta[b]/a) #we vary the dimension of the lattice by considering variation in the temperature, so in beta, since a is fixed
        F[b], F_err[b] = montecarlo_switching(N,n_equil,n_sweeps,n_switching,etha,start,a,delta_x)
        F[b] /= beta[b]
        F[b] += fn.free_energy_zero(beta[b],w0) #this is a constant term related to the choice of the harmonic oscillator as a basis
        F_err[b] /= beta[b]

    np.savetxt(output_path + '/temperature.txt',temperature)
    np.savetxt(output_path + '/free_energy.txt',F)
    np.savetxt(output_path + '/free_energy_err.txt',F_err)
