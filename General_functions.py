import numpy as np
from numba import njit

@njit
def initialize_lattice(n,etha,start):
    if start is True:
        x = np.zeros(n)
        for k in range(n):
            x[k] = -etha
        return x
    elif start is False:
        x = np.random.uniform(-etha,etha,size = n)
        x[0] = x[-2]
        x[-1] = x[1]
        return x

@njit
def calculate_action(x,etha,a):
    S = 0.0
    for i in range(1,x.size):
        S += (1/(4*a))*(np.power((x[i]-x[i-1]),2)) + a*np.power((np.power(x[i],2) - np.power(etha,2)),2)
    return S

@njit
def potential_harmonic(x,w0):
    potential = (1/4.0)*np.power(x*w0,2)
    return potential

@njit
def potential_anharmonic(x, etha):
    potential = np.power(np.power(x,2)-np.power(etha,2),2)
    return potential

@njit
def potential_alpha(x, w0, etha, alpha):
    harmonic = potential_harmonic(x, w0)
    anharmonic = potential_anharmonic(x, etha)
    potential = harmonic + alpha * (anharmonic - harmonic)
    return potential

@njit
def metropolis_switching(x,etha,a,delta_x,alpha):
    w0 = 4*etha
    for i in range(1,x.size):
        old_action = (1/(4*a))*(np.power((x[i]-x[i-1]),2)+np.power((x[i+1]-x[i]),2))+a*potential_alpha(x[i],w0,etha,alpha)
        new_x = x[i]+delta_x*(2*np.random.uniform(0.0,1.0)-1.0)
        new_action = (1/(4*a))*(np.power((new_x-x[i-1]),2)+np.power((x[i+1]-new_x),2))+a*potential_alpha(new_x,w0,etha,alpha)
        delta_action_exp = np.exp(old_action-new_action)
        if delta_action_exp > np.random.uniform(0.0,1.0):
            x[i] = new_x
    x[0] = x[-2]
    x[-1] = x[1]
    return x

@njit
def free_energy_zero(Beta,w0):
    f0 = -np.divide(np.log(2.0*np.sinh(Beta*w0/2.0)),Beta)
    return f0

@njit
def average_std(observable,observable2,n_data):
    observable_av = observable/n_data
    observable_err = observable2/np.power(n_data,2)
    observable_err -= np.power(observable_av,2)/n_data
    observable_err = np.sqrt(observable_err)
    return observable_av,observable_err

@njit
def derivative_log(x,x_err,a):
  derivative_log = np.zeros(29)
  derivative_log_err = np.zeros(29)

  for t in range(29):
    derivative_log[t] = -(x[t+1]-x[t])/(x[t]*a)
    derivative_log_err[t] = np.sqrt(np.power(x_err[t+1]/x[t],2)+np.power(x_err[t]*x[t+1]/np.power(x[t],2),2))/a

  return derivative_log, derivative_log_err

@njit
def metropolis(x, a, delta_x, etha):
  for i in range(1,x.size-1):
    S = (1/(4*a))*(np.power((x[i]-x[i-1]),2)+np.power((x[i+1]-x[i]),2)) + a*np.power((np.power(x[i],2) - np.power(etha,2)),2)
    new_x = x[i] + delta_x*(2.0*np.random.uniform(0.0,1.0)-1.0)
    new_S = (1/(4*a))*(np.power((new_x-x[i-1]),2)+np.power((x[i+1]-new_x),2)) + a*np.power((np.power(new_x,2) - np.power(etha,2)),2)
    delta_S = np.exp(S-new_S)
    if(delta_S > np.random.uniform(0.0,1.0)):
        x[i] = new_x
  x[0] = x[-2]
  x[-1] = x[1]
  return x

@njit
def metropolis_cooling(x_cool, a, delta_x, etha):
  for i in range(1,x_cool.size-1):
    S = (1/(4*a))*(np.power((x_cool[i]-x_cool[i-1]),2)+np.power((x_cool[i+1]-x_cool[i]),2)) + a*np.power((np.power(x_cool[i],2) - np.power(etha,2)),2)
    for _ in range(10):
        new_x_cool = x_cool[i] + delta_x*0.1*(2.0*np.random.uniform(0.0,1.0)-1.0)
        new_S = (1/(4*a))*(np.power((new_x_cool-x_cool[i-1]),2)+np.power((x_cool[i+1]-new_x_cool),2)) + a*np.power((np.power(new_x_cool,2) - np.power(etha,2)),2)
        if(new_S < S):
            x_cool[i] = new_x_cool
  x_cool[0] = x_cool[-2]
  x_cool[-1] = x_cool[1]
  return x_cool

def montecarlo_switching(n,n_equil,n_sweeps,n_switching,etha,start,a,delta_x):
    w0 = 5.6
    d_alpha = 1.0/n_switching
    delta_s_alpha = np.zeros(2*n_switching+1)
    delta_s_alpha2 = np.zeros(2*n_switching+1)
    x = initialize_lattice(n,etha,start)
    for i_switching in range(2*n_switching+1):
        if i_switching <= n_switching:
            alpha = i_switching*d_alpha
        else:
            alpha = 2.0 - i_switching*d_alpha
        for _ in range(n_equil):
            metropolis_switching(x,etha,a,delta_x,alpha)
        for _ in range(n_sweeps-n_equil):
            delta_s_alpha_temp = 0.0
            metropolis_switching(x,etha,a,delta_x,alpha)
            potential_diff = potential_anharmonic(x,etha)
            potential_diff -= potential_harmonic(x,w0)
            delta_s_alpha_temp = np.sum(potential_diff)*a
            delta_s_alpha[i_switching] += delta_s_alpha_temp
            delta_s_alpha2[i_switching] += np.power(delta_s_alpha_temp,2)
    delta_s_alpha_av, delta_s_alpha_err = stat_av_var(delta_s_alpha,delta_s_alpha2,n_sweeps-n_equil)
    integral = simps(delta_s_alpha_av[0:n_switching+1],dx=d_alpha)
    integral_inv = simps(delta_s_alpha_av[n_switching:2*n_switching+1],dx=d_alpha)
    integral_err = simps(delta_s_alpha_err[0:n_switching+1],dx=d_alpha)
    integral_inv_err = simps(delta_s_alpha_err[n_switching:2*n_switching+1],dx=d_alpha)
    propagation_error = np.sqrt(integral_err+integral_inv_err)/2.0
    trapezoidal_error = np.abs((delta_s_alpha_av[n_switching]-delta_s_alpha_av[n_switching-1]+delta_s_alpha_av[1]-delta_s_alpha_av[0])/(d_alpha*12*n_switching*n_switching))
    hysteresis_error = np.abs(integral-integral_inv)
    result = -(integral+integral_inv)/2.0
    result_err = np.sqrt(np.power(propagation_error,2)+np.power(trapezoidal_error,2)+np.power(hysteresis_error,2))
    return result, result_err

@njit
def find_instantons(x,a):
    n_i = n_ai = x_i = x_ai = 0
    pos_i = []
    pos_ai = []
    x0 = x[0]
    for h in range(1,x.size-1):
        if np.abs(x[h])<1e-7:
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
