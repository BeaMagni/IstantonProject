import numpy as np
import General_functions as fn

#monte carlo simulation using the cooling method: computation of correlation functions and their logarithmic derivatives and euclidean path
#the parameters are left as an input from keyboard

def main():
    N = int(input("Insert N, the dimension of the lattice: "))
    etha = float(input("Insert the value for etha, the shift of the potential: "))
    w0 = 4*etha
    a = float(input("Insert the value for a, the lattice spacing: "))
    delta_x = 0.5
    n_equil = int(input("Insert the number of Monte Carlo equilibration sweeps: "))
    n_sweeps = int(input("Insert the number of Monte Carlo sweeps (must be bigger than the previous value): "))
    n_cooling_sweeps = int(input("Insert the number of cooling sweeps: "))
    count = 0
    x_cor_sum = np.zeros(30)
    x2_cor_sum = np.zeros(30)
    x3_cor_sum = np.zeros(30)
    x_cor_sum2 = np.zeros(30)
    x2_cor_sum2 = np.zeros(30)
    x3_cor_sum2 = np.zeros(30)
    x_cor_av = np.zeros(30)
    x2_cor_av = np.zeros(30)
    x3_cor_av = np.zeros(30)
    x_cor_err = np.zeros(30)
    x2_cor_err = np.zeros(30)
    x3_cor_err = np.zeros(30)
    der_log_x  = np.zeros(29)
    der_log_x_err = np.zeros(29)
    der_log_x2  = np.zeros(29)
    der_log_x2_err = np.zeros(29)
    der_log_x3  = np.zeros(29)
    der_log_x3_err = np.zeros(29)

    start = bool(input("Insert the desired start (0-cold, 1-hot): "))
    output_path = './instanton_project/cooling'
    fn.path_creation(output_path)
    
    x = fn.initialize_lattice(N, etha, start)
    
    for _ in range(n_equil): #equilibration runs for the metropolis algorithm
        x = fn.metropolis(x, a, delta_x, etha)
    
    for j in range(n_sweeps-n_equil):
        x = fn.metropolis(x, a, delta_x, etha)
        if j%int((n_sweeps-n_equil)/2) == 0: #we only save some of the configurations
            np.savetxt(output_path + "/montecarlo.txt", x)
        if j%10 == 0: #we don't do the cooling for all the metropolis sweeps
            x_cool = np.copy(x)
            for _ in range(int(n_cooling_sweeps/2)): #we save the configurations only considering 100 cooling sweeps
                x_cool = fn.metropolis_cooling(x_cool, a, delta_x, etha)
            if j%int((n_sweeps-n_equil)/2) == 0: #we only save some of the configurations
                np.savetxt(output_path + "/cooledmontecarlo.txt",x_cool)
            for _ in range(int(n_cooling_sweeps/2),n_cooling_sweeps): #for the correlators we consider 200 cooling sweeps so we continue the cycle
                x_cool = fn.metropolis_cooling(x_cool, a, delta_x, etha)
            for _ in range(5):
                count += 1
                p0 = int((N-30)*np.random.uniform(0.0,1.0)) #we consider a rescaled (based on the total number of lattice points) random number
                for p in range(30):
                    x_cor = x_cool[p0]*x_cool[p0+p]
                    x_cor_sum[p] += x_cor
                    x2_cor_sum[p] += np.power(x_cor,2)
                    x3_cor_sum[p] += np.power(x_cor,3)
                    x_cor_sum2[p] += np.power(x_cor,2)
                    x2_cor_sum2[p] += np.power(x_cor,4)
                    x3_cor_sum2[p] += np.power(x_cor,6)
    
    x1_cor_av, x1_cor_err = fn.average_std(x_cor_sum, x_cor_sum2, count)
    x2_cor_av, x2_cor_err = fn.average_std(x2_cor_sum, x2_cor_sum2, count)
    x3_cor_av, x3_cor_err = fn.average_std(x3_cor_sum, x3_cor_sum2, count)
    der_log_x1, der_log_x1_err = fn.derivative_log(x1_cor_av, x1_cor_err, a)
    x2_cor_av_sub = x2_cor_av-x2_cor_av[-1] #subtraction of the constant term
    x2_cor_err_sub = np.sqrt(np.power(x2_cor_err,2)+np.power(x2_cor_err[-1],2))
    der_log_x2, der_log_x2_err = fn.derivative_log(x2_cor_av_sub, x2_cor_err, a)
    der_log_x3, der_log_x3_err = fn.derivative_log(x3_cor_av, x3_cor_err, a)

    np.savetxt(output_path + '/x1_cor_av.txt',x1_cor_av)
    np.savetxt(output_path + '/x2_cor_av.txt',x2_cor_av)
    np.savetxt(output_path + '/x3_cor_av.txt',x3_cor_av)
    np.savetxt(output_path + '/x1_cor_err.txt',x1_cor_err)
    np.savetxt(output_path + '/x2_cor_err.txt',x2_cor_err)
    np.savetxt(output_path + '/x3_cor_err.txt',x3_cor_err)
    np.savetxt(output_path + '/der_log_x1.txt',der_log_x1)
    np.savetxt(output_path + '/der_log_x1_err.txt',der_log_x1_err)
    np.savetxt(output_path + '/der_log_x2.txt',der_log_x2)
    np.savetxt(output_path + '/der_log_x2_err.txt',der_log_x2_err)
    np.savetxt(output_path + '/der_log_x3.txt',der_log_x3)
    np.savetxt(output_path + '/der_log_x3_err.txt',der_log_x3_err)
