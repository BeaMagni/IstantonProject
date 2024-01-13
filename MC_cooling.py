import numpy as np
from numba import njit
import General_functions as fn

N = 800
etha = 1.4
w0 = 5.6
a = 0.05
delta_x = 0.5
n_equil = 100
n_sweeps = 100000
n_cooling_sweeps = 100
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

x = fn.initialize_lattice(N,etha,start)

for _ in range(n_equil):
    x = fn.metropolis(x,a,delta_x,etha)
for j in range(n_sweeps-n_equil):
    x = fn.metropolis(x,a,delta_x,etha)
    if j%int((n_sweeps-n_equil)/2) == 0:
        np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/montecarlo.txt", x)
    if j%10 == 0:
        x_cool = np.copy(x)
        for _ in range(n_cooling_sweeps):
            x_cool = fn.metropolis_cooling(x_cool,a,delta_x,etha)
        if j%int((n_sweeps-n_equil)/2) == 0:
            np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/cooledmontecarlo.txt",x_cool)
        for _ in range(5):
            count += 1
            p0 = int((N-30)*np.random.uniform(0.0,1.0))
            for p in range(30):
                x_cor = x_cool[p0]*x_cool[p0+p]
                x_cor_sum[p] += x_cor
                x2_cor_sum[p] += np.power(x_cor,2)
                x3_cor_sum[p] += np.power(x_cor,3)
                x_cor_sum2[p] += np.power(x_cor,2)
                x2_cor_sum2[p] += np.power(x_cor,4)
                x3_cor_sum2[p] += np.power(x_cor,6)

x_cor_av, x_cor_err = fn.average_std(x_cor_sum,x_cor_sum2,count)
x2_cor_av, x2_cor_err = fn.average_std(x2_cor_sum,x2_cor_sum2,count)
x3_cor_av, x3_cor_err = fn.average_std(x3_cor_sum,x3_cor_sum2,count)
der_log_x, der_log_x_err = fn.derivative_log(x_cor_av,x_cor_err,a)
x2_cor_av_sub = x2_cor_av-x2_cor_av[-1]
x2_cor_err_sub = np.sqrt(np.power(x2_cor_err,2)+np.power(x2_cor_err[-1],2))
der_log_x2, der_log_x2_err = fn.derivative_log(x2_cor_av_sub,x2_cor_err,a)
der_log_x3, der_log_x3_err = fn.derivative_log(x3_cor_av,x3_cor_err,a)

np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x_cor_av.txt",x_cor_av)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x2_cor_av.txt",x2_cor_av)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x3_cor_av.txt",x3_cor_av)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x_cor_err.txt",x_cor_err)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x2_cor_err.txt",x2_cor_err)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x3_cor_err.txt",x3_cor_err)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x.txt",der_log_x)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x_err.txt",der_log_x_err)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x2.txt",der_log_x2)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x2_err.txt",der_log_x2_err)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x3.txt",der_log_x3)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x3_err.txt",der_log_x3_err)
