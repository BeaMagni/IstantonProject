import numpy as np
from numba import njit #we import this function to speed up the executions of the functions defined in this file
from pathlib import Path

#this function creates the lattice, depending on the choice of hot or cold start

@njit
def initialize_lattice(n, etha, start):
    if start is True:
        x = np.zeros(n)
        for k in range(n):
            x[k] = -etha
        return x
    elif start is False:
        x = np.random.uniform(-etha,etha,size = n) #we generate a random value, uniformly distributed
        x[0] = x[-2] #we impose periodic boundary conditions
        x[-1] = x[1]
        return x

#this function evaluates the discretized action for the anharmonic potential

@njit
def calculate_action(x, etha, a):
    S = 0.0
    for i in range(1,x.size):
        S += (1/(4*a))*(np.power((x[i]-x[i-1]),2)) + a*np.power((np.power(x[i],2) - np.power(etha,2)),2)
    return S

#evaluation of the harmonic potential

@njit
def potential_harmonic(x, w0):
    potential = (1/4.0)*np.power(x*w0,2)
    return potential

#evaluation of the anharmonic potential

@njit
def potential_anharmonic(x, etha):
    potential = np.power(np.power(x,2)-np.power(etha,2),2)
    return potential

#evaluation of the potential used in the adiabatic switching calculations

@njit
def potential_alpha(x, w0, etha, alpha):
    harmonic = potential_harmonic(x, w0)
    anharmonic = potential_anharmonic(x, etha)
    potential = harmonic+alpha*(anharmonic-harmonic)
    return potential

#definition of the metropolis algorithm for the adiabatic switching method

@njit
def metropolis_switching(x, etha, a, delta_x, alpha):
    w0 = 4*etha
    for i in range(1,x.size):
        old_action = (1/(4*a))*(np.power((x[i]-x[i-1]),2)+np.power((x[i+1]-x[i]),2))+a*potential_alpha(x[i],w0,etha,alpha)
        new_x = x[i]+delta_x*(2*np.random.uniform(0.0,1.0)-1.0) #generation of an update for the coordinates, based on a random value
        new_action = (1/(4*a))*(np.power((new_x-x[i-1]),2)+np.power((x[i+1]-new_x),2))+a*potential_alpha(new_x,w0,etha,alpha)
        delta_action_exp = np.exp(old_action-new_action)
        if delta_action_exp > np.random.uniform(0.0,1.0): #acceptance condition
            x[i] = new_x
    x[0] = x[-2] #we impose periodic boundary conditions
    x[-1] = x[1]
    return x

#definition of the constant term to be added to the free energy values

@njit
def free_energy_zero(Beta, w0):
    f0 = -np.divide(np.log(2.0*np.sinh(Beta*w0/2.0)),Beta)
    return f0

#calculations of the statistical averages and errors 

@njit
def average_std(observable, observable2, n_data):
    observable_av = observable/n_data
    observable_err = observable2/np.power(n_data,2)
    observable_err -= np.power(observable_av,2)/n_data
    observable_err = np.sqrt(observable_err)
    return observable_av,observable_err

#computation of the logarithmic derivative of the correlation functions and the corresponding error

@njit
def derivative_log(x, x_err, a):
    derivative_log = np.zeros(29)
    derivative_log_err = np.zeros(29)

    for t in range(29):
        derivative_log[t] = -(x[t+1]-x[t])/(x[t]*a)
        derivative_log_err[t] = np.sqrt((np.power(x_err[t+1]/x[t],2)+np.power(x_err[t]*x[t+1]/np.power(x[t],2),2))/a)

    return derivative_log, derivative_log_err

##definition of the metropolis algorithm

@njit
def metropolis(x, a, delta_x, etha):
    for i in range(1,x.size-1):
        S = (1/(4*a))*(np.power((x[i]-x[i-1]),2)+np.power((x[i+1]-x[i]),2)) + a*np.power((np.power(x[i],2) - np.power(etha,2)),2)
        new_x = x[i] + delta_x*(2.0*np.random.uniform(0.0,1.0)-1.0) #generation of an update for the coordinates, based on a random value
        new_S = (1/(4*a))*(np.power((new_x-x[i-1]),2)+np.power((x[i+1]-new_x),2)) + a*np.power((np.power(new_x,2) - np.power(etha,2)),2)
        delta_S = np.exp(S-new_S)
        if(delta_S > np.random.uniform(0.0,1.0)): #acceptance condition
            x[i] = new_x
    x[0] = x[-2] #we impose periodic boundary conditions
    x[-1] = x[1]
    return x

#definition of the metropolis algorithm for the cooling case

@njit
def metropolis_cooling(x_cool, a, delta_x, etha):
    for i in range(1,x_cool.size-1):
        S = (1/(4*a))*(np.power((x_cool[i]-x_cool[i-1]),2)+np.power((x_cool[i+1]-x_cool[i]),2)) + a*np.power((np.power(x_cool[i],2) - np.power(etha,2)),2)
        for _ in range(10):
            new_x_cool = x_cool[i] + delta_x*0.1*(2.0*np.random.uniform(0.0,1.0)-1.0) #generation of an update for the coordinates, based on a random value
            new_S = (1/(4*a))*(np.power((new_x_cool-x_cool[i-1]),2)+np.power((x_cool[i+1]-new_x_cool),2)) + a*np.power((np.power(new_x_cool,2) - np.power(etha,2)),2)
            if(new_S < S): #acceptance condition: we accept ony value that decrease the action
                x_cool[i] = new_x_cool
    x_cool[0] = x_cool[-2]
    x_cool[-1] = x_cool[1]
    return x_cool

#this functions is used to obtain the instanton - anti-instanton content of the system

@njit
def find_instantons(x, a):
    n_i = n_ai = 0 #instanton and anti-instanton number
    x_i = x_ai = 0
    pos_i = [] #instanton position
    pos_ai = [] #anti-instanton position
    x0 = x[0]
    for h in range(1,x.size-1):
        if np.abs(x[h])<1e-9: #we verify if the value is close to zero
            if x0 > 0:
                n_ai += 1
                x_ai = -(x0*a/(x[h]-x0))+a*(h-1)
                pos_ai.append(x_ai)
            elif x0 < 0:
                n_i += 1
                x_i = -(x0*a/(x[h]-x0))+a*(h-1)
                pos_i.append(x_i)
            else:
                continue
        elif x0*x[h] < 0:
            if x0 > x[h]:
                n_ai += 1
                x_ai = -(x0*a/(x[h]-x0))+a*(h-1)
                pos_ai.append(x_ai)
            elif x0 < x[h]:
                n_i += 1
                x_i = -(x0*a/(x[h]-x0))+a*(h-1)
                pos_i.append(x_i)
        x0 = x[h]
        
    return n_i,n_ai,pos_i,pos_ai

#creation of a directory, in the case it doesn't exist

def path_creation(output_path):
    path = Path(output_path)
    if not path.is_dir():
        path.mkdir(parents=True, exist_ok=True)
