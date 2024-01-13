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
n_cooling_sweeps = 10

start = bool(input("Insert the desired start (0-cold, 1-hot): "))

file_hist = open("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/histogram.txt", 'w')
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

