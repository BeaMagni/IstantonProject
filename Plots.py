import matplotlib.pyplot as plt
import numpy as np

#PLOT POTENTIAL AND ENERGY LEVELS

def potential_eigenvalues():
    etha = 1.4
    x = np.linspace(-2.5,2.5,1000)
    V = (x**2-etha**2)**2
    plt.ylim(0,10)
    plt.xlabel("x")
    plt.ylabel("V(x)")
    plt.plot(x, V, '-')
    energy_eigenvalues = np.loadtxt("./instanton_project/exact/eigenvalues.txt")
    for n in range(4):
        plt.axhline(y=energy_eigenvalues[n], color='red', linestyle='--')
    plt.show()

#PLOT ENERGY EIGENVALUES VARIATION WITH RESPECT TO ETHA

def energy_eigenvalues_variation():
    Etha = np.linspace(0.001,2,1000)
    E_1 = np.loadtxt("./instanton_project/exact/energy_1.txt")
    E_2 = np.loadtxt("./instanton_project/exact/energy_2.txt")
    E_3 = np.loadtxt("./instanton_project/exact/energy_3.txt")
    E_4 = np.loadtxt("./instanton_project/exact/energy_4.txt")
    E_5 = np.loadtxt("./instanton_project/exact/energy_5.txt")
    E_6 = np.loadtxt("./instanton_project/exact/energy_6.txt")
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
    
#PLOT GROUND STATE WAVEFUNCTION HISTOGRAM
def ground_state_hist():
    x = np.loadtxt("./instanton_project/monte_carlo/montecarlo_hist.txt")
    X = np.linspace(-2.5,2.5,100)
    ground_state = np.loadtxt("./instanton_project/exact/ground_state.txt")
    
    plt.plot(X, ground_state, color='lime')
    plt.hist(x, bins=100, density=True, histtype='step')
    plt.xlim(-2.1,2.1)
    plt.ylim(0,0.4)
    plt.xlabel('x')
    plt.ylabel(r'$|\Psi(x)|^2$')
    plt.show()

#PLOT CORRELATION FUNCTIONS

def correlation_functions():
    x1 = np.loadtxt("./instanton_project/monte_carlo/x1_cor_av.txt")
    x2 = np.loadtxt("./instanton_project/monte_carlo/x2_cor_av.txt")
    x3 = np.loadtxt("./instanton_project/monte_carlo/x3_cor_av.txt")
    x1_err = np.loadtxt("./instanton_project/monte_carlo/x1_cor_err.txt")
    x2_err = np.loadtxt("./instanton_project/monte_carlo/x2_cor_err.txt")
    x3_err = np.loadtxt("./instanton_project/monte_carlo/x3_cor_err.txt")
    X1 = np.loadtxt("./instanton_project/exact/x1_corr.txt")
    X2 = np.loadtxt("./instanton_project/exact/x2_corr.txt")
    X3 = np.loadtxt("./instanton_project/exact/x3_corr.txt")
    tau = np.linspace(0,1.5,30)
    Tau = np.linspace(0, 2.5, 100)
    
    plt.xlim(0,1.5)
    plt.ylim(0,8)
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$\langle x^n(0) x^n(\tau) \rangle$')
    plt.plot(Tau, X1, 'b')
    plt.plot(Tau, X2, 'r')
    plt.plot(Tau, X3, 'g')
    plt.errorbar(tau, x1, yerr=x1_err, fmt='o', color='blue', markerfacecolor='none', label='n=1')
    plt.errorbar(tau, x2, yerr=x2_err, fmt='o', color='red', markerfacecolor='none', label='n=2')
    plt.errorbar(tau, x3, yerr=x3_err, fmt='o', color='green', markerfacecolor='none', label='n=3')
    plt.legend()
    plt.show()

#PLOT LOG CORRELATION FUNCTIONS

def log_correlation_functions():
    x1 = np.loadtxt("./instanton_project/monte_carlo/der_log_x1.txt")
    x2 = np.loadtxt("./instanton_project/monte_carlo/der_log_x2.txt")
    x3 = np.loadtxt("./instanton_project/monte_carlo/der_log_x3.txt")
    x1_err = np.loadtxt("./instanton_project/monte_carlo/der_log_x1_err.txt")
    x2_err = np.loadtxt("./instanton_project/monte_carlo/der_log_x2_err.txt")
    x3_err = np.loadtxt("./instanton_project/monte_carlo/der_log_x3_err.txt")
    X1 = np.loadtxt("./instanton_project/exact/log_x1_corr.txt")
    X2 = np.loadtxt("./instanton_project/exact/log_x2_corr.txt")
    X3 = np.loadtxt("./instanton_project/exact/log_x3_corr.txt")
    tau = np.linspace(0,1.5,29)
    Tau = np.linspace(0, 2.5, 100)
    
    plt.xlim(0,1.5)
    plt.ylim(0,5)
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$\frac{d}{d\tau} \log\left\langle x^n(0) x^n(\tau) \right\rangle$')
    plt.plot(Tau, X1, 'b')
    plt.plot(Tau, X2, 'r')
    plt.plot(Tau, X3, 'g')
    plt.errorbar(tau, x1, yerr=x1_err, fmt='o', color='blue', markerfacecolor='none', label='n=1')
    plt.errorbar(tau, x2, yerr=x2_err, fmt='o', color='red', markerfacecolor='none', label='n=2')
    plt.errorbar(tau, x3, yerr=x3_err, fmt='o', color='green', markerfacecolor='none', label='n=3')
    plt.legend()
    plt.show()

#PLOT COOLED CORRELATION FUNCTIONS

def correlation_functions_cooling():
    x1 = np.loadtxt("./instanton_project/cooling/x1_cor_av.txt")
    x2 = np.loadtxt("./instanton_project/cooling/x2_cor_av.txt")
    x3 = np.loadtxt("./instanton_project/cooling/x3_cor_av.txt")
    x1_err = np.loadtxt("./instanton_project/cooling/x1_cor_err.txt")
    x2_err = np.loadtxt("./instanton_project/cooling/x2_cor_err.txt")
    x3_err = np.loadtxt("./instanton_project/cooling/x3_cor_err.txt")
    X1 = np.loadtxt("./instanton_project/exact/x1_corr.txt")
    X2 = np.loadtxt("./instanton_project/exact/x2_corr.txt")
    X3 = np.loadtxt("./instanton_project/exact/x3_corr.txt")
    tau = np.linspace(0,1.5,30)
    Tau = np.linspace(0, 2.5, 100)
    
    plt.xlim(0,1.5)
    plt.ylim(0,8)
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$\langle x^n(0) x^n(\tau) \rangle$')
    plt.plot(Tau, X1, 'b')
    plt.plot(Tau, X2, 'r')
    plt.plot(Tau, X3, 'g')
    plt.errorbar(tau, x1, yerr=x1_err, fmt='o', color='blue', markerfacecolor='none', label='n=1')
    plt.errorbar(tau, x2, yerr=x2_err, fmt='o', color='red', markerfacecolor='none', label='n=2')
    plt.errorbar(tau, x3, yerr=x3_err, fmt='o', color='green', markerfacecolor='none', label='n=3')
    plt.legend()
    plt.show()

#PLOT COOLED LOG CORRELATION FUNCTIONS

def log_correlation_functions_cooling():
    x1 = np.loadtxt("./instanton_project/cooling/der_log_x1.txt")
    x2 = np.loadtxt("./instanton_project/cooling/der_log_x2.txt")
    x3 = np.loadtxt("./instanton_project/cooling/der_log_x3.txt")
    x1_err = np.loadtxt("./instanton_project/cooling/der_log_x1_err.txt")
    x2_err = np.loadtxt("./instanton_project/cooling/der_log_x2_err.txt")
    x3_err = np.loadtxt("./instanton_project/cooling/der_log_x3_err.txt")
    X1 = np.loadtxt("./instanton_project/exact/log_x1_corr.txt")
    X2 = np.loadtxt("./instanton_project/exact/log_x2_corr.txt")
    X3 = np.loadtxt("./instanton_project/exact/log_x3_corr.txt")
    tau = np.linspace(0,1.5,29)
    Tau = np.linspace(0, 2.5, 100)
    
    plt.xlim(0,1.5)
    plt.ylim(0,5)
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$\frac{d}{d\tau} \log\left\langle x^n(0) x^n(\tau) \right\rangle$')
    plt.plot(Tau, X1, 'b')
    plt.plot(Tau, X2, 'r')
    plt.plot(Tau, X3, 'g')
    plt.errorbar(tau, x1, yerr=x1_err, fmt='o', color='blue', markerfacecolor='none', label='n=1')
    plt.errorbar(tau, x2, yerr=x2_err, fmt='o', color='red', markerfacecolor='none', label='n=2')
    plt.errorbar(tau, x3, yerr=x3_err, fmt='o', color='green', markerfacecolor='none', label='n=3')
    plt.legend()
    plt.show()

#PLOT TYPICAL EUCLIDEAN PATH

def euclidean_path():
    x1 = np.loadtxt("./instanton_project/cooling/montecarlo.txt")
    x2 = np.loadtxt("./instanton_project/cooling/cooledmontecarlo.txt")
    tau = np.linspace(0, 40, 800)
    
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

def free_energy():
    temperature = np.loadtxt("./instanton_project/adiabatic/temperature.txt")
    Temperature = np.linspace(0.01,2,100)
    Free_energy = np.loadtxt("./instanton_project/exact/free_energy.txt")
    free_energy = np.loadtxt("./instanton_project/adiabatic/free_energy.txt")
    free_energy_err = np.loadtxt("./instanton_project/adiabatic/free_energy_err.txt")
    
    plt.ylim(-2.5,-1.0)
    plt.xlim(0.02,2.5)
    plt.ylabel('F')
    plt.xlabel('T')
    plt.plot(Temperature, Free_energy, color='black', label='Exact')
    plt.errorbar(temperature, free_energy, yerr=free_energy_err, fmt='o', color='blue', markerfacecolor='none', label='Adiabatic switching')
    plt.xscale('log')
    plt.legend()
    plt.show()

#PLOT INSTANTON - ANTI-INSTANTON DISTRIBUTION HISTOGRAM

def instanton_distribution():
    zcr = np.loadtxt("./instanton_project/instantons/histogram.txt", float, delimiter = ' ')
    plt.hist(zcr, 40, (0.0,4.0), histtype='step', color='black')
    plt.xlim(0,4)
    plt.ylim(0,4000)
    plt.xlabel(r'$\tau_z$')
    plt.ylabel(r'$n_{IA}(\tau_z)$')
    plt.show()

#PLOT INSTANTON DENSITY vs NUMBER OF COOLING SWEEPS

def instanton_density():
    n_cooling = np.loadtxt("./instanton_project/instantons/n_cooling.txt")
    
    colors = ['red', 'green', 'blue']
    i = 0
    
    for etha in range(14,17,1):
        etha = etha/10
        n_tot = np.loadtxt(f"./instanton_project/instantons/n_tot_av_{etha}.txt")
        n_tot_err = np.loadtxt(f"./instanton_project/instantons/n_tot_err_{etha}.txt")
        plt.errorbar(n_cooling, n_tot, yerr=n_tot_err, fmt='o', color=colors[i], markerfacecolor='none', label=f'$\eta = {etha}$')
        s0 = (4/3)*np.power(etha,3)
        loop_1 = 8*np.power(etha,5/2)*np.sqrt(2/np.pi)*np.exp(-s0)
        loop_2 = 8*np.power(etha,5/2)*np.sqrt(2/np.pi)*np.exp(-s0-71/(72*s0))
        plt.axhline(y=loop_1, color=colors[i])
        plt.axhline(y=loop_2, color=colors[i], linestyle='--')
        i += 1
    
    plt.xlabel(r'$n_{cool}$')
    plt.ylabel(r'$N_{inst}/\beta$')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(top=1)
    plt.xlim(1,200)
    plt.legend()
    plt.show()

#PLOT INSTANTON ACTION DENSITY vs NUMBER OF COOLING SWEEPS

def instanton_action():
    n_cooling = np.loadtxt("./instanton_project/instantons/n_cooling.txt")
    
    colors = ['red', 'green', 'blue']
    i = 0
    
    for etha in range(14,17,1):
        etha = etha/10
        action = np.loadtxt(f"./instanton_project/instantons/action_av_{etha}.txt")
        action_err = np.loadtxt(f"./instanton_project/instantons/action_err_{etha}.txt")
        plt.errorbar(n_cooling, action, yerr=action_err, fmt='o', color=colors[i], markerfacecolor='none', label=f'$\eta = {etha}$')
        s0 = (4/3)*np.power(etha,3)
        plt.axhline(y=s0, color=colors[i])
        i += 1
    
    plt.xlabel(r'$n_{cool}$')
    plt.ylabel(r'$S/N_{inst}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1,200)
    plt.legend()
    plt.show()
