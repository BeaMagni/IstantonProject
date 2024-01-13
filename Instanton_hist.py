import numpy as np
from numba import njit
import General_functions as fn

def main():
    N = int(input("Insert N, the dimension of the lattice: "))
    etha = float(input("Insert the value for etha, the shift of the potential: "))
    w0 = 4*etha
    a = float(input("Insert the value for a, the lattice spacing: "))
    delta_x = 0.5
    n_equil = int(input("Insert the number of Monte Carlo equilibration sweeps: "))
    n_sweeps = int(input("Insert the number of Monte Carlo sweeps (must be bigger than the previous value): "))
    n_cooling_sweeps = int(input("Insert the number of cooling sweeps: "))
    output _path = './instanton_project/instantons'
    fn.path_creation(output_path)
    start = bool(input("Insert the desired start (0-cold, 1-hot): "))
    
    file_hist = open(output_path + '/histogram.txt', 'w')
    x = fn.initialize_lattice(N,etha,start)
    
    for _ in range(n_equil):
        x = fn.metropolis(x,a,delta_x,etha)
    for j in range(n_sweeps-n_equil):
        x = fn.metropolis(x,a,delta_x,etha)
        if j%10 == 0:
            x_cool = np.copy(x)
            n_inst = n_anti_inst = 0
            for _ in range(n_cooling_sweeps):
                x_cool = fn.metropolis_cooling(x_cool,a,delta_x,etha)
            n_inst, n_anti_inst, inst_pos, anti_inst_pos = fn.find_instantons(x_cool,a)
            zia = 0
            if n_inst >0 and n_inst == n_anti_inst:
                if inst_pos[0] < anti_inst_pos[0]:
                    for k in range(n_inst):
                        if k == 0:
                            zm = anti_inst_pos[-1] - N*a
                        else:
                            zm = anti_inst_pos[k-1]
                        zia = np.minimum(np.abs(anti_inst_pos[k]-inst_pos[k]),np.abs(inst_pos[k]-zm))
                        file_hist.write(str(zia)+'\n')
                elif inst_pos[0] > anti_inst_pos[0]:
                    for k in range(n_inst):
                        if k == 0:
                            zp = inst_pos[-1] - N*a
                        else:
                            zp = inst_pos[k-1]
                        zia = np.minimum(np.abs(-anti_inst_pos[k]+inst_pos[k]),np.abs(anti_inst_pos[k]-zp))
                        file_hist.write(str(zia)+'\n')
                else:
                    continue

