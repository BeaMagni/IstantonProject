
from asyncore import read
import numpy as np
import math
from numpy import linalg as la
import matplotlib.pyplot as plt
from scipy.special import factorial
from numpy.polynomial import hermite


#building of the hamiltonian of the anharmonic oscillator and computation of energy eigenvalues E and eigenvector v
def diag(ndim, etha):
    
    H = np.zeros((ndim,ndim))

    for n in range(ndim):
        x = 3 * A * (c**4) * ((n+1)**2+n**2)+B * (c**2) * (2*n+1)+om_0 * (n+1/2)+C
        H[n][n] = x
        if n+2<ndim:
          x_2 = A*c**4*(4*n+6)*np.sqrt((n+1)*(n+2))+B*c**2*np.sqrt((n+1)*(n+2))
          H[n][n+2] = H[n+2][n] = x_2
        if n+4<ndim:
          x_4 = c**4*np.sqrt((n+1)*(n+2)*(n+3)*(n+4))
          H[n][n+4] = H[n+4][n] = x_4

    E, v = la.eigh(H)
    return E,v

#behaviour of the energy eigenvalues varying etha 
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
      om_0 = 4*etha
      B = -(2*etha**2+(om_0**2)/4)
      C = etha**4
      c = 1/np.sqrt(om_0)
      H = np.zeros((ndim,ndim))
      for n in range(ndim):
        x = 3 * A * (c**4) * ((n+1)**2+n**2)+B * (c**2) * (2*n+1)+om_0 * (n+1/2)+C
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

    plt.ylim(0,22)
    plt.xlim(0.25,2)
    plt.xlabel(r'$\eta$')
    plt.ylabel("E")
    plt.yticks(np.arange(0,25,5))
    plt.plot(Etha, vector1, 'b-')
    plt.plot(Etha, vector2, 'b-')
    plt.plot(Etha, vector3, 'b-')
    plt.plot(Etha, vector4, 'b-')
    plt.plot(Etha, vector5, 'b-')
    plt.plot(Etha, vector6, 'b-')
    return plt.show()


def en_densities(ndim, etha, v):
    om_0= 4*etha
    c = 1/np.sqrt(om_0)
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
    return corr_funct

def log_correlation(ndim, etha, tau, ntau, tau_array, dtau, E, v):
    log_corr_funct = np.zeros((3,ntau))
    corr_funct= correlation_fuctions(ndim, etha, tau, ntau, tau_array, dtau, E, v)
    i = j = 0
    energy_densities = en_densities(ndim, etha, v)
    for tau in tau_array:
      for rho in energy_densities[0]:
        log_corr_funct[0][i] += (rho*(E[j]-E[0])*np.exp(-(E[j]-E[0])*tau))/(corr_funct[0][i])
        j += 1
      j = 0
      i += 1

    i = j = 0

    for tau in tau_array:
      for rho in energy_densities[1]:
        log_corr_funct[1][i] += (rho*(E[j]-E[0])*np.exp(-(E[j]-E[0])*tau))/(corr_funct[1][i]-energy_densities[1,0])
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
    return log_corr_funct

def free_energy(E):
    temperature = np.linspace(0.01, 2.0, 999)

    free_energy = []
    Z= 0.0
    i_f_energy = 0
    for t in temperature:
        for e in E:
            Z += np.exp(-e / t)

        free_energy.append(t * np.log(Z))
        Z = 0.0
    plt.yticks([-2.5,-2,-1.5,-1])
    plt.xlim(0.05,2.5)
    plt.ylim(-2.5,-1)
    plt.xscale('log')
    plt.ylabel('F')
    plt.xlabel('T')
    plt.plot(temperature, free_energy)
    plt.show()

def hermite_coefficients(position,norm,ndim):
  coefficient = np.zeros(ndim)
  for n in range(ndim):
    coefficient[n] = pow(np.sqrt(np.pi)*norm*math.factorial(n)*pow(2.0,n),-1/2)*np.exp(-0.5*pow(position/norm,2))
  return coefficient

def psi_ground_state(position, norm, ndim, v):
    ground_state = np.zeros(position.size)
    for x in range(position.size):
        projections = np.multiply(v[:,0],hermite_coefficients(X[x],norm,ndim))
        ground_state[x] = pow(hermite.hermval(position[x]/norm,projections),2)
    plt.xlim(-2.1,2.1)
    plt.ylim(0,0.4)
    plt.xlabel('x')
    plt.ylabel(r'$|\Psi(x)|^2$')
    plt.plot(position,ground_state,'b-')
    plt.show()

#main
ndim = int(input("Insert the dimention of the matrix: "))
etha = float(input("Insert the value for the shift of the potential: "))
om_0 = c = C = B = 0
A = 1
om_0 = 4*etha
B = -2*etha**2-om_0**2/4
C = etha**4
c = 1/np.sqrt(om_0)
energy_eigenvalues, energy_eigenvectors = diag(ndim,etha)
x = np.linspace(-2.5,2.5,1000)
V=(x**2-etha**2)**2

plt.ylim(0,10)
plt.xlabel("x")
plt.ylabel("V(x)")
plt.plot(x, V, '-')


for n in range(4):
  plt.axhline(y=energy_eigenvalues[n], color='red', linestyle='--')

plt.show()

energy_variation(ndim)

euclidian_time = 2.5
step = 100
t_array, dt = np.linspace(0, euclidian_time, step, retstep=True)
corr = correlation_fuctions(ndim, etha, euclidian_time, step, t_array, dt, energy_eigenvalues, energy_eigenvectors)
log_corr= log_correlation(ndim, etha, euclidian_time, step, t_array, dt, energy_eigenvalues, energy_eigenvectors)
plt.xlim(0,1.5)
plt.ylim(0,8)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\langle x^n(0) x^n(\tau) \rangle$')
plt.plot(t_array,corr[0],'r', label='n=1')
plt.plot(t_array,corr[1],'b', label='n=2')
plt.plot(t_array,corr[2],'g', label='n=3')
plt.legend()
plt.show()
plt.xlim(0,1.5)
plt.ylim(0,5)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\frac{d}{d\tau} \log\left\langle x^n(0) x^n(\tau) \right\rangle$')
plt.plot(t_array,log_corr[0], 'r', label='n=1')
plt.plot(t_array,log_corr[1], 'b', label='n=2')
plt.plot(t_array,log_corr[2], 'g', label='n=3')
plt.legend()
plt.show()
free_energy(energy_eigenvalues)

X = np.linspace(-2.5,2.5,100)
norm = c*np.sqrt(2)
psi_ground_state(X,norm,ndim, energy_eigenvectors)