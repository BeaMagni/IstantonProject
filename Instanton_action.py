import numpy as np
from numba import njit
import FUNCTIONS as fn

N = 800
a = 0.05
delta_x = 0.5
n_equil = 100
n_sweeps = 100000
n_cooling_sweeps = 200
count = n_cooling = n_tot = 0

start = bool(input("Insert the desired start (0-cold, 1-hot): "))

for etha in range(14,17,1):
    etha = etha/10
    x = fn.initialize_lattice(N,etha,start)
    n_tot_sum = np.zeros(n_cooling_sweeps)
    n2_tot_sum = np.zeros(n_cooling_sweeps)
    action_cooling = np.zeros(n_cooling_sweeps)
    action2_cooling = np.zeros(n_cooling_sweeps)
    n_cooling = 0
    for _ in range(n_equil):
        x = fn.metropolis(x,a,delta_x,etha)
    for j in range(n_sweeps-n_equil):
        x = fn.metropolis(x,a,delta_x,etha)
        if j%10 == 0:
            x_cool = np.copy(x)
            n_cooling += 1
            n_inst = n_anti_inst = 0
            for u in range(n_cooling_sweeps):
                x_cool = fn.metropolis_cooling(x_cool,a,delta_x,etha)
                n_inst, n_anti_inst, _, _ = fn.find_instantons(x_cool,a)
                n_tot = n_inst+n_anti_inst
                n_tot_sum[u] += n_tot
                n2_tot_sum[u] += np.power(n_tot,2)
                action = fn.calculate_action(x_cool,etha,a)
                action_cooling[u] += action
                action2_cooling[u] += np.power(action,2)
    action_av, action_err = fn.stat_av_var(action_cooling, action2_cooling,n_cooling)
    n_tot_av, n_tot_err = fn.stat_av_var(n_tot_sum,n2_tot_sum,n_cooling)
    action_av = np.divide(action_av,n_tot_av)
    action_err = np.divide(action_err,n_tot_av)
    n_tot_av /= N*a
    n_tot_err /= N*a
    np.savetxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/n_tot_av_{etha}.txt",n_tot_av)
    np.savetxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/n_tot_err_{etha}.txt",n_tot_err)
    np.savetxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/action_av_{etha}.txt",action_av)
    np.savetxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/action_err_{etha}.txt",action_err)

x_axis = np.linspace(1,n_cooling_sweeps+1,n_cooling_sweeps,False)
np.savetxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/n_cooling.txt", x_axis)
