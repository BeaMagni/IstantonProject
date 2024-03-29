import numpy as np
from numba import njit
import General_functions as fn

#monte carlo simulation: computation of correlation functions and their logarithmic derivatives and the groundstate wavefunction
#the parameters are left as an input from keyboard

def main():
    N = int(input("Insert N, the dimension of the lattice: "))
    etha = float(input("Insert the value for etha, the shift of the potential: "))
    w0 = 4*etha
    a = float(input("Insert the value for a, the lattice spacing: "))
    delta_x = 0.5
    n_equil = int(input("Insert the number of Monte Carlo equilibration sweeps: "))
    n_sweeps = int(input("Insert the number of Monte Carlo sweeps (must be bigger than the previous value): "))
    x_cor_sum = np.zeros(30)
    x2_cor_sum = np.zeros(30)
    x3_cor_sum = np.zeros(30)
    x_cor_sum2 = np.zeros(30)
    x2_cor_sum2 = np.zeros(30)
    x3_cor_sum2 = np.zeros(30)
    x_cor_av = np.zeros(30)
    x2_cor_av = np.zeros(30)
    x3_cor_av = np.zeros(30)
    x_cor2_av = np.zeros(30)
    x2_cor2_av = np.zeros(30)
    x3_cor2_av = np.zeros(30)
    x_cor_err = np.zeros(30)
    x2_cor_err = np.zeros(30)
    x3_cor_err = np.zeros(30)
    der_log_x  = np.zeros(29)
    der_log_x_err = np.zeros(29)
    der_log_x2  = np.zeros(29)
    der_log_x2_err = np.zeros(29)
    der_log_x3  = np.zeros(29)
    der_log_x3_err = np.zeros(29)
    count = 0
    
    start = bool(input("Insert the desired start (0-cold, 1-hot): "))
    output_path = './instanton_project/monte_carlo'
    fn.path_creation(output_path)
    
    x = fn.lattice_initialization(N, etha, start)
    
    for _ in range(n_equil): #equilibration runs for the metropolis algorithm
        x = fn.metropolis(x, a, delta_x, etha)
    with open(output_path + "/montecarlo_hist.txt", 'w') as hist:    
        for j in range(n_sweeps - n_equil):
          x = fn.metropolis(x, a, delta_x, etha)
          if j%50 == 0: #we only save some of the configurations to avoid correlation between various configurations
            np.savetxt(hist, x)
          for _ in range(5):
            count += 1
            p0 = int((N-30)*np.random.uniform(0.0,1.0)) #we consider a rescaled (based on the total number of lattice points) random number
            for p in range(30):
              x_cor = x[p0]*x[p0+p]
              x_cor_sum[p] += x_cor
              x2_cor_sum[p] += np.power(x_cor,2)
              x3_cor_sum[p] += np.power(x_cor,3)
              x_cor_sum2[p] += np.power(x_cor,2)
              x2_cor_sum2[p] += np.power(x_cor,4)
              x3_cor_sum2[p] += np.power(x_cor,6)
        
    x_cor_av, x_cor_err = fn.average_std(x_cor_sum, x_cor_sum2, count)
    x2_cor_av, x2_cor_err = fn.average_std(x2_cor_sum, x2_cor_sum2, count)
    x3_cor_av, x3_cor_err = fn.average_std(x3_cor_sum, x3_cor_sum2, count)
    der_log_x, der_log_x_err = fn.derivative_log(x_cor_av, x_cor_err, a)
    x2_cor_av_sub = x2_cor_av-x2_cor_av[-1] #subtraction of the constant term
    x2_cor_err_sub = np.sqrt(np.power(x2_cor_err,2)+np.power(x2_cor_err[-1],2))
    der_log_x2, der_log_x2_err = fn.derivative_log(x2_cor_av_sub, x2_cor_err_sub, a)
    der_log_x3, der_log_x3_err = fn.derivative_log(x3_cor_av, x3_cor_err, a)

    np.savetxt(output_path + '/x1_cor_av.txt',x_cor_av)
    np.savetxt(output_path + '/x2_cor_av.txt',x2_cor_av)
    np.savetxt(output_path + '/x3_cor_av.txt',x3_cor_av)
    np.savetxt(output_path + '/x1_cor_err.txt',x_cor_err)
    np.savetxt(output_path + '/x2_cor_err.txt',x2_cor_err)
    np.savetxt(output_path + '/x3_cor_err.txt',x3_cor_err)
    np.savetxt(output_path + '/der_log_x1.txt',der_log_x)
    np.savetxt(output_path + '/der_log_x1_err.txt',der_log_x_err)
    np.savetxt(output_path + '/der_log_x2.txt',der_log_x2)
    np.savetxt(output_path + '/der_log_x2_err.txt',der_log_x2_err)
    np.savetxt(output_path + '/der_log_x3.txt',der_log_x3)
    np.savetxt(output_path + '/der_log_x3_err.txt',der_log_x3_err)
    
