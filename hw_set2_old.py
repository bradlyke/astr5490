import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
matplotlib.rc('text',usetex=True)
import sys

def table_load():
    #This will load the table of the orbital parameters
    orbit_tab = np.loadtxt('orbit_params.dat',skiprows=3,dtype=bytes,delimiter=',').astype(str)
    col_names = np.array(['NAME','a','e','i','w','O','P'])
    angle_names = np.array(['i','w','O'])
    col_forms = np.array(['U7','f4','f4','f4','f4','f4','f4'])
    orbit_arr = np.zeros(8,dtype={'names':col_names,'formats':col_forms})
    orbit_arr['NAME'] = orbit_tab[:,0]
    for i in range(1,7):
        orbit_arr[col_names[i]] = orbit_tab[:,i].astype('float32')
    for cname in angle_names:
        orbit_arr[cname] = orbit_arr[cname] * (np.pi/180.0) #converts degrees to radians
    orbit_arr['P'] = orbit_arr['P'] * 3.15576e7 #converts period to seconds
    return orbit_arr

#This is part of the function for part c. It will find E_i+1 from E_i
def e_next(E_p,M,e):
    denom = (1 - (e*np.cos(E_p)))**(-1.0) #derivative of g(E)
    numer = E_p - (e*np.sin(E_p)) - M #g(E)
    E_next = E_p - (numer * denom) #The definition of Newton-Raphson method
    return E_next

#This is the second part for part c. It iterates through E until
#E_i+1 - E_i is less than the check_val.
#It runs through a number of steps equal to counter to keep a diverging solution
#from running forever. If it fails it returns a fail_check of 1
def newt_method(M,e,E0):
    counter = 10000 #Number of steps for protection
    fail_check = 0
    check_val = 1e-5 #Threshhold for convergence. Set smaller for better solution
    E_p = E0 #Initialize E_i for the first step
    for i in range(counter):
        E_n = e_next(E_p,M,e)
        if np.abs(E_n - E_p) <= check_val: #Check for final solution.
            E_final = E_n
            break
        else: #If it hasn't converged yet, go to the next iteration.
            E_p = E_n
        if i +1 == counter: #This is the check in case it diverged by the step limit.
            print('Did not converge in {} steps'.format(counter))
            fail_check = 1
            E_final = E_n
    #Return some value for E_final (hopefully the convergent one) and the
    #fail_check condition.
    return E_final,fail_check

#This is first function for part d. Finds radius (r) from E.
def r_func(E_f,e,a):
    return a*(1 - (e*np.cos(E_f)))

#This is the second function for part d. Finds the angle from pericenter
#based on r_current and the input eccentricity and semi-major axis.
def f_func(r,e,a):
    numer = (a * (1-e**(2))) - r
    denom = (e * r)**(-1.0)
    f = np.arccos(numer * denom) #Uses multiplication to keep from blowing up
    return f

def pos_vel(a,e,i,wom,Om,T,t):
    #This will find the instantaneous velocity and position
    pos_vec = np.zeros(3,dtype='f8') #0 - X, 1 - Y, 2 - Z
    vel_vec = np.zeros(3,dtype='f8') #0 - Vx, 1 - Vy, 2 - Vz
    #Calculate r (radius) first
    n = (2*np.pi) / T
    M = n * t
    E_instant,fc = newt_method(M,e,0)
    r = r_func(E_instant,e,a)
    f = f_func(r,e,a)
    xang = ((np.cos(Om)*np.cos(wom+f)) - (np.sin(Om)*np.sin(wom+f)*np.cos(i)))
    yang = ((np.sin(Om)*np.cos(wom+f)) - (np.cos(Om)*np.sin(wom+f)*np.cos(i)))
    zang = np.sin(wom+f)*np.sin(i)
    pos_vec[0] = r*xang
    pos_vec[1] = r*yang
    pos_vec[2] = r*zang

    rdot = (a*e*(1-e**(2))*np.sin(f)) * (1 + (e*np.cos(f)))**(-1.0)
    vel_vec[0] = (rdot*xang) - (r*((np.cos(Om)*np.sin(wom+f)) + (np.sin(Om)*np.cos(wom+f)*np.cos(i))))
    vel_vec[1] = (rdot*yang) - (r*((np.sin(Om)*np.sin(wom+f)) - (np.cos(Om)*np.cos(wom+f)*np.cos(i))))
    vel_vec[2] = (rdot*zang) + (r*np.cos(wom+f)*np.sin(i))

    return pos_vec,vel_vec


if __name__=='__main__':
    tab_test = table_load()
    #test_planet = input('Please type planet name: ')
    test_planet = 'earth'
    wplanet = np.where(tab_test['NAME']==test_planet)[0]
    a0 = tab_test['a'][wplanet][0]
    e0 = tab_test['e'][wplanet][0]
    i0 = tab_test['i'][wplanet][0]
    wom0 = tab_test['w'][wplanet][0]
    Om0 = tab_test['O'][wplanet][0]
    T0 = tab_test['P'][wplanet][0]
    orbit_vec = np.array([a0,e0,i0,wom0,Om0,T0])
    num_steps = 1000
    step_size = T0/num_steps
    pos_array = np.zeros((num_steps,3),dtype='f8')
    vel_array = np.zeros((num_steps,3),dtype='f8')
    for i in range(num_steps):
        t_step = i*step_size
        pos_array[i,:],vel_array[i,:] = pos_vel(*orbit_vec,t_step)
    back_half_check = 
    fig,ax = plt.subplots(figsize=(10,8))
    ax.scatter(pos_array[:,0],pos_array[:,1],color='black',marker='.')
    plt.show()
