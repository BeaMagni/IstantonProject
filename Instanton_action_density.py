import numpy as np
import General_functions as fn

#computation of the instanton density and instanton action as functions of the number of cooling sweeps. 
#the parameters are left as an input from keyboard

def main():
    N = int(input("Insert N, the dimension of the lattice: "))
    a = float(input("Insert the value for a, the lattice spacing: "))
    delta_x = 0.5
    n_equil = int(input("Insert the number of Monte Carlo equilibration sweeps: "))
    n_sweeps = int(input("Insert the number of Monte Carlo sweeps (must be bigger than the previous value): "))
    n_cooling_sweeps = int(input("Insert the number of cooling sweeps: "))
    count = n_cooling = n_tot = 0
    output_path = './instanton_project/instantons'
    fn.path_creation(output_path)
    start = bool(input("Insert the desired start (0-cold, 1-hot): "))

    #we consider three values of etha to evaluate the instanton desnity and action
    
    for etha in range(14,17,1):
        etha = etha/10 
        x = fn.lattice_inizialization(N,etha,start)
        n_tot_sum = np.zeros(n_cooling_sweeps)
        n2_tot_sum = np.zeros(n_cooling_sweeps)
        action_cooling = np.zeros(n_cooling_sweeps)
        action2_cooling = np.zeros(n_cooling_sweeps)
        n_cooling = 0
        
        for _ in range(n_equil): #equilibration runs for the metropolis algorithm
            x = fn.metropolis(x, a, delta_x, etha)
        
        for j in range(n_sweeps-n_equil):
            x = fn.metropolis(x, a, delta_x, etha)
            if j%10 == 0: #we don't do the cooling for all the metropolis sweeps
                x_cool = np.copy(x)
                n_cooling += 1
                n_inst = n_anti_inst = 0
                for u in range(n_cooling_sweeps):
                    x_cool = fn.metropolis_cooling(x_cool, a, delta_x, etha)
                    n_inst, n_anti_inst, _, _ = fn.instantons_content(x_cool, a) #we only consider the number of instantons and anti-instantons
                    n_tot = n_inst+n_anti_inst
                    n_tot_sum[u] += n_tot
                    n2_tot_sum[u] += np.power(n_tot,2)
                    action = fn.calculate_action(x_cool, etha, a)
                    action_cooling[u] += action
                    action2_cooling[u] += np.power(action,2)
        action_av, action_err = fn.average_std(action_cooling, action2_cooling, n_cooling)
        n_tot_av, n_tot_err = fn.average_std(n_tot_sum, n2_tot_sum, n_cooling)
        action_av = np.divide(action_av,n_tot_av)
        action_err = np.divide(action_err,n_tot_av)
        n_tot_av /= N*a
        n_tot_err /= N*a
        np.savetxt(output_path + f'/n_tot_av_{etha}.txt',n_tot_av)
        np.savetxt(output_path + f'/n_tot_err_{etha}.txt',n_tot_err)
        np.savetxt(output_path + f'/action_av_{etha}.txt',action_av)
        np.savetxt(output_path + f'/action_err_{etha}.txt',action_err)
    
    x_axis = np.linspace(1,n_cooling_sweeps+1,n_cooling_sweeps,False)
    np.savetxt(output_path + '/n_cooling.txt', x_axis)
