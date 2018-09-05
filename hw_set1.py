"""
By: Brad Lyke
All of the code for homework set 1. This requires python 3 and latex (for plots).

Parameters
----------
None : :class:'float'
	The program will prompt for user input.
    On the active menu a character is need.
    For each subprogram inputs must be floats.

Returns
----------
--hello--
:class:'str'
	hello world string

--lin_interp--
:class:'float'
    Returns the y value for a input x value (based on two points).

--newt_method--
:class:'float','int'
    Returns the convergent value for E and a fail_check condition.

--orbit_plot--
:class:'int'
    Returns a fail_check condition. If it passes, returns two plots.
"""
import numpy as np
import sys
import os

#This function will format an answer string in a box to stand out.
def ans_box(ans_str):
    box_str = '#'*(len(ans_str) + 4) #Make the box string the right length
    print('\n')
    print(box_str)
    print('# {} #'.format(ans_str)) #Print the actual answer
    print(box_str)

#This is the program for part a.
def hello():
    h_str = 'hello world'
    return h_str

#This is the function for part b. It assumes the two points define a line.
#Takes two points and the text X value.
def lin_interp(x1,y1,x2,y2,x_t):
    m = (y2 - y1) / (x2 - x1) #Find the slope
    y = m*(x_t - x1) + y1 #Use the slope to find y_test with x_test
    return y

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

#This is the final function for part d. Actually makes the plots of the orbital
#path and f as a function of time. It requires three user inputs, the mean motion
#in rad/sec, the eccentricity, and the semi-major axis.
#For clearer plots, set the number of steps (num_steps) higher.
#Returns the plot for f(t) for t in days, but calculates in seconds.
#Requires: 0<e<1
#If newton's method above fails, then this function will fail and print a
#statement of such.
def orbit_plot(n,e,a):
    import matplotlib.pyplot as plt
    import matplotlib
    import matplotlib.ticker as ticker
    matplotlib.rc('text',usetex=True)

    num_steps = 100
    T_sec = (2*np.pi) / n #calculate the period in seconds from mean motion
    step_size = T_sec / num_steps #number of seconds per step
    #Set up an array of the parameters as they are calculated. Pre-allocates
    #the array on the number of steps for optimization.
    param_array = np.zeros(num_steps,dtype=[('t','f4'),('M','f4'),
                            ('E','f4'),('r','f4'),('f','f4')])
    #The time per step and mean anomaly can be predefined. Assume that t0 = 0
    param_array['t'] = np.array([i*step_size for i in range(num_steps)])
    param_array['M'] = n*param_array['t']
    #Assume initial position is at periastron. This means a phase = 0.
    param_array['r'][0] = a * (1-e)
    fail_check2 = 0 #fail_check for this function.
    #iterate through time steps to find E,r, and f.
    for i in range(1,num_steps):
        param_array['E'][i],fail_check = newt_method(param_array['M'][i],e,param_array['E'][i-1])
        #If newton's method for eccentric anomaly fails, orbital path fails.
        if fail_check == 1:
            fail_check2 = 1
            print('E diverged, bad values')
            break
        #If it didn't fail, find r(t) and f(t).
        param_array['r'][i] = r_func(param_array['E'][i],e,a)
        param_array['f'][i] = f_func(param_array['r'][i],e,a)
        #If f(t) is exactly pi, then arccos fails and returns NaN. This fixes
        #that edge case.
        if np.isnan(param_array['f'][i]):
            param_array['f'][i] = np.pi
        #Arccos only returns values between (-1,1]. This mirrors around x-axis
        #so that the orbital path works correctly.
        back_half_check = param_array['t'][i] / T_sec
        if (back_half_check > 0.5):
            param_array['f'][i] += 2*(np.pi - param_array['f'][i])
    #If we mapped the whole orbital path, then plot the path and plot f(t).
    if fail_check2 == 0:
        x_arr = param_array['r'] * np.cos(param_array['f'])
        y_arr = param_array['r'] * np.sin(param_array['f'])
        fig1,ax1 = plt.subplots(figsize=(5,4))
        ax1.scatter(x_arr,y_arr,color='black',marker='.')
        ax1.set_xlabel('X-Coordinate (m)',fontsize=15)
        ax1.set_ylabel('Y-Coordinate (m)',fontsize=15)
        ax1.tick_params(axis='both',direction='in')
        ax1.tick_params(axis='both',which='minor',direction='in')
        ax1.tick_params(top=True,right=True)
        ax1.tick_params(which='minor',top=True,right=True)
        #f(t) will be plotted with t in days (seconds is ridiculous)
        fig2,ax2 = plt.subplots(figsize=(5,4))
        ax2.scatter(param_array['t']/86400,param_array['f'],color='black',marker='.')
        ax2.set_xlabel(r'$t$ (days)',fontsize=15)
        ax2.set_ylabel(r'$f$ (rad)',fontsize=15)
        ax2.tick_params(axis='both',direction='in')
        ax2.tick_params(axis='both',which='minor',direction='in')
        ax2.tick_params(top=True,right=True)
        ax2.tick_params(which='minor',top=True,right=True)
        plt.show()
        return 0
    elif fail_check2 == 1:
        return 1

#A user might want to run more than one function in a row or run a function more
#than once. This is the active menu that will do that.
#All subfunctions other than hello world take multiple user inputs. These NEED
#to be comma separated.
def active_menu():
    while True:
        print("\nWelcome to Brad's Homework Set 1")
        print('--------------------------------\n')
        print('Please select a subprogram')
        print('---NOTE: separate multiple inputs with , and no spaces')
        print('[a] Hello World\n[b] Linear Interpolation')
        print('[c] Root Finding\n[d] Orbital Plot\n[e] Quit')
        user_choice = input('Selection: ')
        if ((user_choice == 'a')or(user_choice == '1')):
            h_ans = hello()
            ans_box(h_ans)
            print('\n')
        elif ((user_choice == 'b')or(user_choice == '2')):
            print('\nThis program will do a linear interpolation between two points')
            x_a_str,y_a_str = input('Input first point (x_1,y_1): ').split(',')
            x_a,y_a = float(x_a_str),float(y_a_str)
            x_b_str,y_b_str = input('Input second point (x_2,y_2): ').split(',')
            x_b,y_b = float(x_b_str),float(y_b_str)
            x_test = float(input('Input x-coordinate to test: '))
            y_test = lin_interp(x_a,y_a,x_b,y_b,x_test)
            y_ans_str = 'Y_test = {:.2f}'.format(y_test)
            ans_box(y_ans_str)
            print('\n')
        elif ((user_choice == 'c')or(user_choice == '3')):
            print("\nThis program will solve Kepler's Eqn. for E")
            print('Input M and E_init below in rad/s and rad respectively')
            M_s,e_s,E_init = input('Please input values (M,e,E_init): ').split(',')
            M,e,E0 = float(M_s),float(e_s),float(E_init)
            E_f,fail_check = newt_method(M,e,E0)
            if fail_check == 0:
                g_check = E_f - (e*np.sin(E_f)) - M
                E_fstr = 'E_final = {}'.format(E_f)
                g_str = 'g_check = {}'.format(g_check)
                ans_box(E_fstr)
                print("\ng_check won't be exactly 0, but this should be close:")
                ans_box(g_str)
            else:
                print('Bad Values')
                print('\n')
        elif ((user_choice == 'd')or(user_choice == '4')):
            print('\nThis program will solve for the orbital path of an exoplanet')
            print('Input n in rad/s and a in meters')
            n_str,e_str,a_str = input('Please input values (n,e,a): ').split(',')
            n,e,a = float(n_str),float(e_str),float(a_str)
            fc2 = orbit_plot(n,e,a)
            if fc2 == 1:
                print('Try running with different values.')
            print('\n')
        else:
            print('\nBye!')
            break

#Will run the active menu if the program is run from the command line.
if __name__ == '__main__':
    active_menu()
