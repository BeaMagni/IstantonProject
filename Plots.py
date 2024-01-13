import matplotlib.pyplot as plt
import numpy as np

#PLOT POTENSTIAL AND ENERGY LEVELS
def potential_eigenvalues():
    x = np.linspace(-2.5,2.5,1000)
    V=(x**2-etha**2)**2
    plt.ylim(0,10)
    plt.xlabel("x")
    plt.ylabel("V(x)")
    plt.plot(x, V, '-')
    energy_eigenvalues = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/energy_eigenvalues.txt")
    for n in range(4):
        plt.axhline(y=energy_eigenvalues[n], color='red', linestyle='--')
    plt.show()

#PLOT ENERGY VARIATION
def energy__eigenvalues_variation():
    Etha = np.linspace(0.001,2,1000)
    E_1 = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/energy_1.txt")
    E_2 = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/energy_2.txt")
    E_3 = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/energy_3.txt")
    E_4 = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/energy_4.txt")
    E_5 = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/energy_5.txt")
    E_6 = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/energy_6.txt")
    plt.ylim(0,22)
    plt.xlim(0.25,2)
    plt.xlabel(r'$\eta$')
    plt.ylabel("E")
    plt.yticks(np.arange(0,25,5))
    plt.plot(Etha, E_1, 'b-')
    plt.plot(Etha, E_2, 'b-')
    plt.plot(Etha, E_3, 'b-')
    plt.plot(Etha, E_4, 'b-')
    plt.plot(Etha, E_5, 'b-')
    plt.plot(Etha, E_6, 'b-')
    plt.show()
    
#PLOT GROUND STATE HISTOGRAM
def ground_state_hist():
    x = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/hist.txt")
    y = np.linspace(-2.5,2.5,100)
    ground_state = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/groundstate.txt")
    
    plt.plot(y,ground_state,color='lime',label='Exact')
    plt.hist(x, bins=100, density=True, histtype='step',label='Monte Carlo')
    plt.xlim(-2.1,2.1)
    plt.ylim(0,0.4)
    plt.xlabel('x')
    plt.ylabel(r'$|\Psi(x)|^2$')
    plt.show()


#PLOT CORRELATION FUNCTIONS

x = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/x_cor_av.txt")
x2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/x2_cor_av.txt")
x3 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/x3_cor_av.txt")
x_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/x_cor_std.txt")
x2_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/x2_cor_std.txt")
x3_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/x3_cor_std.txt")
X = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/x_cor.txt")
X2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/x2_cor.txt")
X3 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/x3_cor.txt")
tau = np.linspace(0,1.5,30)
Tau = np.linspace(0, 2.5, 100)

plt.xlim(0,1.5)
plt.ylim(0,8)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\langle x^n(0) x^n(\tau) \rangle$')
plt.plot(Tau,X, 'b')
plt.plot(Tau,X2, 'r')
plt.plot(Tau,X3, 'g')
plt.errorbar(tau,x,yerr=x_err,fmt='o', color='blue', markerfacecolor='none', label='n=1')
plt.errorbar(tau,x2,yerr=x2_err,fmt='o', color='red', markerfacecolor='none', label='n=2')
plt.errorbar(tau,x3,yerr=x3_err,fmt='o', color='green', markerfacecolor='none', label='n=3')
plt.legend()
plt.show()


#PLOT LOG CORRELATION FUNCTIONS

x = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/der_log_x.txt")
x2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/der_log_x2.txt")
x3 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/der_log_x3.txt")
x_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/der_log_x_err.txt")
x2_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/der_log_x2_err.txt")
x3_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/montecarlo/der_log_x3_err.txt")
X = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/log_x_cor.txt")
X2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/log_x2_cor.txt")
X3 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/log_x3_cor.txt")
tau = np.linspace(0,1.5,29)
Tau = np.linspace(0, 2.5, 100)

plt.xlim(0,1.5)
plt.ylim(0,5)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\frac{d}{d\tau} \log\left\langle x^n(0) x^n(\tau) \right\rangle$')
plt.plot(Tau,X, 'b')
plt.plot(Tau,X2, 'r')
plt.plot(Tau,X3, 'g')
plt.errorbar(tau,x,yerr=x_err,fmt='o', color='blue', markerfacecolor='none', label='n=1')
plt.errorbar(tau,x2,yerr=x2_err,fmt='o', color='red', markerfacecolor='none', label='n=2')
plt.errorbar(tau,x3,yerr=x3_err,fmt='o', color='green', markerfacecolor='none', label='n=3')
plt.legend()
plt.show()



#PLOT COOLED CORRELATION FUNCTIONS

x = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x_cor_av.txt")
x2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x2_cor_av.txt")
x3 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x3_cor_av.txt")
x_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x_cor_std.txt")
x2_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x2_cor_std.txt")
x3_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/x3_cor_std.txt")
X = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/x_cor.txt")
X2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/x2_cor.txt")
X3 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/x3_cor.txt")
tau = np.linspace(0,1.5,30)
Tau = np.linspace(0, 2.5, 100)

plt.xlim(0,1.5)
plt.ylim(0,8)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\langle x^n(0) x^n(\tau) \rangle$')
plt.plot(Tau,X, 'b')
plt.plot(Tau,X2, 'r')
plt.plot(Tau,X3, 'g')
plt.errorbar(tau,x,yerr=x_err,fmt='o', color='blue', markerfacecolor='none', label='n=1')
plt.errorbar(tau,x2,yerr=x2_err,fmt='o', color='red', markerfacecolor='none', label='n=2')
plt.errorbar(tau,x3,yerr=x3_err,fmt='o', color='green', markerfacecolor='none', label='n=3')
plt.legend()
plt.show()

#PLOT COOLED LOG CORRELATION FUNCTIONS

x = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x.txt")
x2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x2.txt")
x3 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x3.txt")
x_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x_err.txt")
x2_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x2_err.txt")
x3_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/der_log_x3_err.txt")
X = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/log_x_cor.txt")
X2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/log_x2_cor.txt")
X3 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/log_x3_cor.txt")
tau = np.linspace(0,1.5,29)
Tau = np.linspace(0, 2.5, 100)

plt.xlim(0,1.5)
plt.ylim(0,5)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\frac{d}{d\tau} \log\left\langle x^n(0) x^n(\tau) \right\rangle$')
plt.plot(Tau,X, 'b')
plt.plot(Tau,X2, 'r')
plt.plot(Tau,X3, 'g')
plt.errorbar(tau,x,yerr=x_err,fmt='o', color='blue', markerfacecolor='none', label='n=1')
plt.errorbar(tau,x2,yerr=x2_err,fmt='o', color='red', markerfacecolor='none', label='n=2')
plt.errorbar(tau,x3,yerr=x3_err,fmt='o', color='green', markerfacecolor='none', label='n=3')
plt.legend()
plt.show()



#PLOT CONFIGURATIONS

x1 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/montecarlo.txt")
x2 = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/cooling/cooledmontecarlo.txt")
tau = np.linspace(0, 40,800)

plt.plot(tau, x1, color='black', linewidth=0.5, label='Monte Carlo')
plt.plot(tau, x2, color='lime', linewidth=2, label='Cooled Monte Carlo')
plt.yticks([-2,-1,0,1,2])
plt.xticks([0,5,10,15,20])
plt.xlim(0,20)
plt.ylim(-2,2)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$x(\tau)$')
plt.legend()
plt.show()



#PLOT FREE ENERGY

temperature = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/adiabatic/temperature.txt")
Temperature = np.linspace(0.01,2,100)
Free_energy = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/exact/free_energy.txt")
free_energy = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/adiabatic/free_energy.txt")
free_energy_err = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/adiabatic/free_energy_err.txt")

plt.ylim(-2.5,-1.0)
plt.xlim(0.02,2.5)
plt.ylabel('F')
plt.xlabel('T')
plt.plot(Temperature,Free_energy, color='black', label='Exact')
plt.errorbar(temperature,free_energy,yerr=free_energy_err,fmt='o', color='blue', markerfacecolor='none', label='Adiabatic switching')
plt.xscale('log')
plt.legend()
plt.show()



#INSTANTON HISTOGRAM

zcr = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/histogram.txt", float, delimiter = ' ')
plt.hist(zcr,40,(0.0,4.0),histtype='step', color='black')
plt.xlim(0,4)
#plt.ylim(0,4000)
plt.xlabel(r'$\tau_z$')
plt.ylabel(r'$n_{IA}(\tau_z)$')

plt.show()



#INSTANTON DENSITY

n_cooling = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/n_cooling.txt")

colors = ['red', 'green', 'blue']
i = 0

for etha in range(14,17,1):
    etha = etha/10
    n_tot = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/n_tot_av_{etha}.txt")
    n_tot_err = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/n_tot_err_{etha}.txt")
    plt.errorbar(n_cooling,n_tot,n_tot_err,fmt='o', color=colors[i], markerfacecolor='none', label=f'$\eta = {etha}$')
    s0 = (4/3)*np.power(etha,3)
    loop_1 = 8*np.power(etha,5/2)*np.sqrt(2/np.pi)*np.exp(-s0)
    loop_2 = 8*np.power(etha,5/2)*np.sqrt(2/np.pi)*np.exp(-s0-71/(72*s0))
    plt.axhline(y=loop_1,color=colors[i])
    plt.axhline(y=loop_2,color=colors[i],linestyle='--')
    i +=1

plt.xlabel(r'$n_{cool}$')
plt.ylabel(r'$N_{inst}/\beta$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(top=1)
plt.xlim(1,200)
plt.legend()
plt.show()



#INSTANTON ACTION DENSITY

n_cooling = np.loadtxt("C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/n_cooling.txt")

colors = ['red', 'green', 'blue']
i = 0

for etha in range(14,17,1):
    etha = etha/10
    action = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/action_av_{etha}.txt")
    action_err = np.loadtxt(f"C:/Users/115271/Desktop/UniBO/Theoretical and Numerical Aspects of Nuclear Physics/esame/instantons/action_err_{etha}.txt")
    plt.errorbar(n_cooling,action,action_err,fmt='o', color=colors[i], markerfacecolor='none', label=f'$\eta = {etha}$')
    s0 = (4/3)*np.power(etha,3)
    plt.axhline(y=s0,color=colors[i])
    i +=1

plt.xlabel(r'$n_{cool}$')
plt.ylabel(r'$S/N_{inst}$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1,200)
plt.legend()
plt.show()

