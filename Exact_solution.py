import numpy as np
import math
from numpy import linalg as la 
from scipy.special import factorial
from numpy.polynomial import hermite  #this library is needed to compute the eigenfunction for the ground state
import General_functions as fn

output_path = './instanton_project/exact'
fn.path_creation(output_path)

#building of the hamiltonian of the anharmonic oscillator and computation of energy eigenvalues E and eigenvector v

def diag(ndim, etha, A, B, C, c):
    w0 = 4*etha
    H = np.zeros((ndim,ndim))
    f = []

    for n in range(ndim):
        x = 3*A*(c**4)*((n+1)**2+n**2)+B*(c**2)*(2*n+1)+w0*(n+1/2)+C
        H[n][n] = x
        if n+2<ndim:
          x_2 = A*c**4*(4*n+6)*np.sqrt((n+1)*(n+2))+B*c**2*np.sqrt((n+1)*(n+2))
          H[n][n+2] = H[n+2][n] = x_2
        if n+4<ndim:
          x_4 = c**4*np.sqrt((n+1)*(n+2)*(n+3)*(n+4))
          H[n][n+4] = H[n+4][n] = x_4

    E, v = la.eigh(H) #this function diagonalizes the matrix given as an input and returns both eigenvalues (E) and eigenvectors (v)
    for i in range(4):
        f.append(E[i])
    np.savetxt(output_path + '/eigenvalues.txt', f)
    
    return E,v

#evaluation of the energy eigenvalues as the parameter etha varies

def energy_variation(ndim):
    A = 1
    vector1 = []
    vector2 = []
    vector3 = []
    vector4 = []
    vector5 = []
    vector6 = []
    
    Etha = np.linspace(0.001,2,1000)

    for i in range(len(Etha)):
        etha = float(Etha[i])
        w0 = 4*etha
        B = -(2*etha**2+(w0**2)/4)
        C = etha**4
        c = 1/np.sqrt(w0)
        H = np.zeros((ndim,ndim))
        for n in range(ndim):
            x = 3 * A * (c**4) * ((n+1)**2+n**2)+B * (c**2) * (2*n+1)+w0 * (n+1/2)+C
            H[n][n] = x
            if n+2<ndim:
                x_2 = A*c**4*(4*n+6)*np.sqrt((n+1)*(n+2))+B*c**2*np.sqrt((n+1)*(n+2))
                H[n][n+2] = H[n+2][n] = x_2
            if n+4<ndim:
                x_4 = c**4*np.sqrt((n+1)*(n+2)*(n+3)*(n+4))
                H[n][n+4] = H[n+4][n] = x_4
    
        E, _ = la.eigh(H)
        vector1.append(E[0])
        vector2.append(E[1])
        vector3.append(E[2])
        vector4.append(E[3])
        vector5.append(E[4])
        vector6.append(E[5])
        
    np.savetxt(output_path + '/energy_1.txt',vector1)
    np.savetxt(output_path + '/energy_2.txt',vector2)
    np.savetxt(output_path + '/energy_3.txt',vector3)
    np.savetxt(output_path + '/energy_4.txt',vector4)
    np.savetxt(output_path + '/energy_5.txt',vector5)
    np.savetxt(output_path + '/energy_6.txt',vector6)

#evaluation of the energy densities that will be used in computing correlation functions

def en_densities(ndim, etha, v):
    w0= 4*etha
    c = 1/np.sqrt(w0)
    energy_densities = np.zeros((3,ndim))
    for n in range(ndim):
        cn = 0.0
        dn = 0.0
        en = 0.0
        for m in range(ndim):
            mm3 = max(m-3,0)
            mm2 = max(m-2,0)
            mm1 = max(m-1,0)
            mp1 = min(m+1,ndim-1)
            mp2 = min(m+2,ndim-1)
            mp3 = min(m+3,ndim-1)
    
            cn += (pow(m+1,1/2)*v[mp1,0]+pow(m,1/2)*v[mm1,0])*v[m,n]
            dn += (pow(m*(m-1),1/2)*v[mm2,0]+(2*m+1)*v[m,0]+pow((m+1)*(m+2),1/2)*v[mp2,0])*v[m,n]
            en += (pow(m*(m-1)*(m-2),1/2)*v[mm3,0]+3*m*pow(m,1/2)*v[mm1,0]+3*(m+1)*pow(m+1,1/2)*v[mp1,0]+pow((m+1)*(m+2)*(m+3),1/2)*v[mp3,0])*v[m,n]
        
        energy_densities[0,n] = pow(c,2)*pow(cn,2)
        energy_densities[1,n] = pow(c,4)*pow(dn,2)
        energy_densities[2,n] = pow(c,6)*pow(en,2)
    
    return energy_densities

#computation of the exact correlation functions for x, x^2 and x^3

def correlation_fuctions(ndim, etha, tau, ntau, tau_array, dtau, E, v):
    corr_funct = np.zeros((3,ntau))
    energy_densities = en_densities(ndim, etha, v)
    i = j = 0

    for tau in tau_array:
        for rho in energy_densities[0]:
            corr_funct[0][i] += rho*np.exp(-(E[j]-E[0])*tau)
            j += 1
        j = 0
        i += 1

    i = j = 0

    for tau in tau_array:
        for rho in energy_densities[1]:
            corr_funct[1][i] += rho*np.exp(-(E[j]-E[0])*tau)
            j += 1
        j = 0
        i += 1

    i = j = 0

    for tau in tau_array:
        for rho in energy_densities[2]:
            corr_funct[2][i] += rho*np.exp(-(E[j]-E[0])*tau)
            j += 1
        j = 0
        i += 1
    
    for i in range(3):
        np.savetxt(output_path + f'/x{i+1}_corr.txt', corr_funct[i]) 
    
    return corr_funct

#computation of the logarithmic derivative of the correlation functions for x, x^2 and x^3

def log_correlation(ndim, etha, tau, ntau, tau_array, dtau, E, v):
    log_corr_funct = np.zeros((3,ntau))
    corr_funct = correlation_fuctions(ndim, etha, tau, ntau, tau_array, dtau, E, v)
    energy_densities = en_densities(ndim, etha, v)
    
    i = j = 0
    for tau in tau_array:
        for rho in energy_densities[0]:
            log_corr_funct[0][i] += (rho*(E[j]-E[0])*np.exp(-(E[j]-E[0])*tau))/(corr_funct[0][i])
            j += 1
        j = 0
        i += 1

    i = j = 0

    for tau in tau_array:
        for rho in energy_densities[1]:
            log_corr_funct[1][i] += (rho*(E[j]-E[0])*np.exp(-(E[j]-E[0])*tau))/(corr_funct[1][i]-energy_densities[1,0]) #for the x^2 correlation function we must remove the constant contribution
            j += 1
        j = 0
        i += 1

    i = j = 0

    for tau in tau_array:
        for rho in energy_densities[2]:
            log_corr_funct[2][i] += (rho*(E[j]-E[0])*np.exp(-(E[j]-E[0])*tau))/(corr_funct[2][i])
            j += 1
        j = 0
        i += 1
        
    for i in range(3):
        np.savetxt(output_path + f"/log_x{i+1}_corr.txt",log_corr_funct[i])

#evaluation of the free energy of the system 

def free_energy(E):
    temperature = np.linspace(0.01, 2.0, 999)
    free_energy = []
    Z= 0.0
    i_f_energy = 0
    for t in temperature:
        for e in E:
            Z += np.exp(-e/t)

        free_energy.append(t*np.log(Z))
        Z = 0.0
    np.savetxt(output_path + '/free_energy.txt', free_energy)

#definition of the coefficients to use in the Hermite polinomials expansion of the ground state wavefunction

def hermite_coefficients(position, norm, ndim):
    coefficient = np.zeros(ndim)
    for n in range(ndim):
        coefficient[n] = pow(np.sqrt(np.pi)*norm*math.factorial(n)*pow(2.0,n),-1/2)*np.exp(-0.5*pow(position/norm,2))
    return coefficient

#computation of the ground state wavefunction

def psi_ground_state(position, norm, ndim, v):
    ground_state = np.zeros(position.size)
    for x in range(position.size):
        projections = np.multiply(v[:,0],hermite_coefficients(X[x],norm,ndim))
        ground_state[x] = pow(hermite.hermval(position[x]/norm,projections),2)
    np.savetxt(output_path + '/ground_state.txt', ground_state)

#main: all the inputs are left as an input from keyboard

def main():
    ndim = int(input("Insert the dimention of the matrix: "))
    etha = float(input("Insert the value for the shift of the potential: "))
    euclidian_time = float(input("Insert the maximum value for the euclidean time: "))
    step = int(input("Insert the value for the time step: "))
    w0 = c = C = B = 0
    A = 1
    w0 = 4*etha
    B = -2*etha**2-w0**2/4
    C = etha**4
    c = 1/np.sqrt(w0)
    t_array, dt = np.linspace(0, euclidian_time, step, retstep=True)
    X = np.linspace(-2.5,2.5,100)
    norm = c*np.sqrt(2)
    
    energy_eigenvalues, energy_eigenvectors = diag(ndim, etha, A, B, C, c)
    energy_variation(ndim)
    
    correlation_fuctions(ndim, etha, euclidian_time, step, t_array, dt, energy_eigenvalues, energy_eigenvectors)
    log_correlation(ndim, etha, euclidian_time, step, t_array, dt, energy_eigenvalues, energy_eigenvectors)
    free_energy(energy_eigenvalues)
    psi_ground_state(X, norm, ndim, energy_eigenvectors)

