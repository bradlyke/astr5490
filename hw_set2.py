"""
By: Brad Lyke
All of the code for homework set 2. This requires python 3 and latex (for plots).
This also requires that sciCon be in the same folder (or in your PYTHONPATH).
It is included in the utilities repository (or submitted along with this program).

Parameters
----------
None : :class:'float'
	The program will prompt for user input.
    For user choices, use 1/2 or y/n as asked.
    For the initial position/velocity input comma-separated (no spaces)
      values in floats or scientific notation (centimeters for R, km/s for V)

Returns
----------
--main--
:class:'str'
    Returns the orbital parameters generated from either the test point or the
    user's inputs.

--prob2_plotter--
:class:'none'
    Returns 3 figures:
        a) (Y v X) and (y v x)
        b) X(t) and Y(t)
        c) Vz(t)

--prob3_plotter--
:class:'none'
    Returns 3 figures:
        a) Five plots for variable values of i (e~0)
        b) Five plots for variable w (e=0.5)
        c) Five plots for variable O (e=0.5,w=pi/2)

Note: The initial position and velocity input by the user must match a
self-consistent set of data or the program may return garbage (which is not
this programmer's responsibility to correct for). Make sure V^(2) will work
for a given R^(2). Also, the programs assume the R,V given are for periastron.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
matplotlib.rc('text',usetex=True)
import sys
from sciCon import mks #I wrote this, so you'll have to have it in the same folder.
#from mpl_toolkits.mplot3d import Axes3D

#This is part of the function for root finding from homework 1.
#It will find E_i+1 from E_i
def e_next(E_p,M,e):
    denom = (1 - (e*np.cos(E_p)))**(-1.0) #derivative of g(E)
    numer = E_p - (e*np.sin(E_p)) - M #g(E)
    E_next = E_p - (numer * denom) #The definition of Newton-Raphson method
    return E_next

#This is the second part for root finding from homework 1.
#It iterates through E until E_i+1 - E_i is less than the check_val.
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

#This is radial function for finding the position. Finds radius (r) from E.
#Note, this and the root functions require e to not be 0. Later code where e=0
#will actually use e=10^-8
def r_func(E_f,e,a):
    return a*(1 - (e*np.cos(E_f)))

#This is the angle function for finding the position.
#Finds the angle from pericenter based on r_current and the
#input eccentricity and semi-major axis.
#Because this function uses arccos, which has a limited range, code calling
#this function will need to correct f(t) if t > T/2.
def f_func(r,e,a):
    numer = (a * (1-e**(2))) - r
    denom = (e * r)**(-1.0)
    f = np.arccos(numer * denom) #Uses multiplication to keep from blowing up
    return f

#This function finds a vector's magnitude. I needed this enough to write a function
def vmag(in_vec):
    return np.sqrt(in_vec.dot(in_vec))

#This function finds the angle from the cross product definition of
#AxB = ABsin(x). Will find the angle in both radians and degrees by default.
def cross_angle(vector1,vector2,deg=True):
    cprod = np.cross(vector1,vector2)
    prod_mag = vmag(cprod)
    vec1_mag = vmag(vector1)
    vec2_mag = vmag(vector2)
    denom = vec1_mag * vec2_mag
    cang = np.arcsin(prod_mag * denom**(-1.0))
    cang_deg = cang * (180/np.pi)
    if deg==True:
        return cang,cang_deg
    else:
        return cang

#This does the same as above, but using the dot product definition
#A.B = ABcos(x)
def dot_angle(vector1,vector2,deg=True):
    dprod = vector1.dot(vector2)
    vec1_mag = vmag(vector1)
    vec2_mag = vmag(vector2)
    denom = vec1_mag * vec2_mag
    dang = np.arccos(dprod * denom**(-1.0))
    dang_deg = dang * (180/np.pi)
    if deg==True:
        return dang,dang_deg
    else:
        return dang

#This solves for the orbital parameters (a, e, Omega, i, w, Period, |h|)
#for a given initial position vector and velocity vector. Because a single point
#is not enough to solve for these, the second point is assumed to be the same
#point in the observer frame, but at periastron. If the point isn't periastron
#too bad.
def orbit_params(obs_init,vel_obs_init):
    #This will find the instantaneous velocity and position
    #Calculate semi-major axis first
    G = mks.G
    mSun = 1.98855e30 #These are in mks. Units don't really matter as long as they
    mJup = 1.8982e27  #are consistent.
    mT = mSun+mJup
    GMT = G*mT
    r_init = np.sqrt(obs_init[0]**(2) + obs_init[1]**(2)+obs_init[2]**(2))
    v_initsq = vel_obs_init[0]**(2)+vel_obs_init[1]**(2)+vel_obs_init[2]**(2)
    a = ((2/r_init) - (v_initsq/(GMT)))**(-1.0)
    #We are assuming we are at periastron, with the initial position and velocity
    pos_init = np.array([r_init,0,0])
    vel_init = np.array([0,np.sqrt(v_initsq),0])
    #Now find the period
    T = np.sqrt((4*np.pi**(2)*a**(3))/(GMT))
    #Now find the magnitude of the angular momentum (per unit mass) in the Orbital
    #frame of reference.
    h = np.cross(pos_init,vel_init)
    hmag = vmag(h)
    #Finally find the eccentricity
    e_part = (T*hmag)/(2*np.pi*a**(2))
    e = np.sqrt(1 - e_part**(2))
    r = a*(1-e)
    #To find the angles for the orbit, we need the angular momentum per unit mass
    #in the observer frame as well, which is H.
    H = np.cross(obs_init,vel_obs_init)
    Hmag = vmag(H)
    Zhat = np.array([0,0,1]) #Observer frame Zhat vector.
    nvec = np.cross(Zhat,H)
    #Find Omega based on H.
    if nvec[1] >= 0:
        Omega = (np.arccos(nvec[0]/vmag(nvec)))
    elif nvec[1] < 0:
        Omega = (2*np.pi) - (np.arccos(nvec[0]/vmag(nvec)))
    #The inclination can be found from Hxh as the only rotation between these two
    #is i, due to h having no z-component.
    inc = cross_angle(H,h,deg=False)
    #The distance from the star in the observer's frame.
    Rmag = vmag(obs_init)
    #Finally, solve the function for Z(t) (observer frame) for womega (the last)
    #unknown angle).
    womega = np.arcsin(obs_init[2]/(Rmag*np.sin(inc)))
    #Print out the orbital parameters we just found and return them.
    print('Semi-Major Axis: {:5.3f} au'.format(a/mks.au))
    print('Eccentricity: {:6.4f}'.format(e))
    print('Longitude of Ascending Node (O): {:6.2f} deg'.format(Omega*180/np.pi))
    print('Inclination (i): {:6.2f} deg'.format(inc*180/np.pi))
    print('Argument of Periapse (w): {:6.2f} deg'.format(womega*180/np.pi))

    return a,e,Omega,inc,womega,T,hmag

#This function will find the instantaneous position and velocity (observer's frame)
#from the orbital parameters and t (time). To ensure the function for E converges
#we need to keep the previous E found from the last iteration and use it for our
#initial E in this iteration.
def pos_vel(a,e,Om,i,wom,t,T,bhc,hmag,Ep):
    #Initialize empty vectors. orb_vec is the position in the orbital plane. Used
    #for plotting purposes in problem 2 (for comparison).
    orb_vec = np.zeros(3,dtype='f8') #0 - x, 1 - y, 2 - z
    pos_vec = np.zeros(3,dtype='f8') #0 - X, 1 - Y, 2 - Z
    vel_vec = np.zeros(3,dtype='f8') #0 - Vx, 1 - Vy, 2 - Vz
    #Find the mean motion.
    n = (2*np.pi) / T
    #Find the mean anomaly based on the time.
    M = n * t
    #Solve for E, then r, then f for a given time.
    E_instant,fc = newt_method(M,e,Ep)
    r = r_func(E_instant,e,a)
    f = f_func(r,e,a)
    #This corrects the angle, f(t), based on the period due to f(t) using arccos.
    if bhc == 1:
        f += 2*(np.pi - f)

    #This finds the orbital plane position. Note that z is always 0.
    orb_vec[0] = r*np.cos(f)
    orb_vec[1] = r*np.sin(f)
    orb_vec[2] = 0

    #Now find the position and velocity in the observer's frame. The angular
    #component of the position reappears in the velocity (from the product rule)
    #so I find it separately. Note that fdot is based on |h|, not |H|, so we need
    #the magnitude of the angular momentum vector in the orbital plane as |h| is
    #constant.
    xang = (np.cos(Om)*np.cos(wom + f)) - (np.sin(Om)*np.sin(wom + f)*np.cos(i))
    yang = (np.sin(Om)*np.cos(wom + f)) + (np.cos(Om)*np.sin(wom + f)*np.cos(i))
    zang = np.sin(wom + f)*np.sin(i)
    pos_vec[0] = r*xang
    pos_vec[1] = r*yang
    pos_vec[2] = r*zang

    #And now we generate the velocity for the given point. Note that rdot depends
    #on fdot, but fdot only relies on |h| and r (which we already found).
    #I took these derivatives by hand.
    fdot = hmag * r**(-2.0)
    rdot = a*e*(1-e**(2))*np.sin(f)*fdot*(1 + (e*np.cos(f)))**(-2.0)
    vel_vec[0] = (rdot*xang) - (r*fdot*((np.cos(Om)*np.sin(wom+f)) + (np.sin(Om)*np.cos(wom+f)*np.cos(i))))
    vel_vec[1] = (rdot*yang) - (r*fdot*((np.sin(Om)*np.sin(wom+f)) - (np.cos(Om)*np.cos(wom+f)*np.cos(i))))
    vel_vec[2] = (rdot*zang) + (r*fdot*np.cos(wom+f)*np.sin(i))

    #Returns the position and velocity for the observer frame. Also returns the
    #position for the orbital plane (for plotting) and the value we found for E.
    return pos_vec,vel_vec,orb_vec,E_instant

#Orbital parameters must be self-consistent. If the speed is too low, the object
#won't orbit. If the distance is too long, the speed might be too high.
#This function generates a self-consistent test point (with a given set of
#parameters) so I could test the above code and ensure it was giving back the
#right values for the angles. If you choose the first option from the program
#this is the function called.
def gen_test_point(a_test,e_test,i_test,w_test,o_test):
    #Convert the values for a, e, i, w, O into the right units (centimeters and radians)
    a = a_test*mks.au
    e = e_test
    i = i_test*np.pi/180
    w = w_test*np.pi/180
    O = o_test*np.pi/180
    #We assume this object is at periastron. If it's not, too bad.
    f = 0
    r = a*(1-e)
    #Find v^2, observer position, and observer velocity vector from givens.
    mSun = 1.98855e30
    mJup = 1.8982e27
    mT = mSun+mJup
    pos_orb_test = np.array([r,0,0]) #periastron orbital-frame position vector.
    vsq = mks.G*mT*((2/r) - (1/a))
    vel_orb_test = np.array([0,np.sqrt(vsq),0]) #Periastron orbital-frame v-vector.
    #Observer-frame position.
    xang = (np.cos(O)*np.cos(w))-(np.sin(O)*np.sin(w)*np.cos(i))
    yang = (np.sin(O)*np.cos(w))+(np.cos(O)*np.sin(w)*np.cos(i))
    zang = np.sin(w)*np.sin(i)
    X0 = r*xang
    Y0 = r*yang
    Z0 = r*zang
    #Note that the functions for V(t), don't include the rdot term Because
    #we don't need an rdot to make it work for this single point, just the rotation
    #into the correct frame.
    VX0 = -np.sqrt(vsq)*((np.cos(O)*np.sin(w)) + (np.sin(O)*np.cos(w)*np.cos(i)))
    VY0 = -np.sqrt(vsq)*((np.sin(O)*np.sin(w)) - (np.cos(O)*np.cos(w)*np.cos(i)))
    VZ0 = np.sqrt(vsq)*np.cos(w)*np.sin(i)
    #Turn the above information into vectors.
    pos_obs_test = np.array([X0,Y0,Z0])
    vel_obs_test = np.array([VX0,VY0,VZ0])

    #Return the vectors found. r0, R0, v0, V0
    return pos_orb_test,pos_obs_test,vel_orb_test,vel_obs_test

#This function creates the plots asked for in problem 2. The first plot is a 3D
#plot of the observer-frame elliptical orbit in X,Y,Z. This requires a non-standard
#package, so I commented it out. If you want to see it, remove the quotes and uncomment
#the line at the top importing Axes3D.
def prob2_plotter(pos_array,orb_array,vel_array,t_arr):
    matplotlib.rc('font',size=15)
    '''
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111,projection='3d')
    ax1.scatter(pos_array[:,0],pos_array[:,1],pos_array[:,2])
    ax1.set_xlabel('X-Coordinate (m)')
    ax1.set_ylabel('Y-Coordinate (m)')
    ax1.set_zlabel('Z-Coordinate (m)')
    ax1.tick_params(axis='both',direction='in')
    ax1.tick_params(axis='both',which='minor',direction='in')
    ax1.tick_params(top=True,right=True)
    ax1.tick_params(which='minor',top=True,right=True)
    '''

    #This plots x,y and X,Y on the same plot so you can see what the rotation
    #does to the ellipse.
    fig2,ax2 = plt.subplots(figsize=(10,8))
    ax2.scatter(orb_array[:,0]/mks.au,orb_array[:,1]/mks.au,color='magenta',marker='.',label='Orbital Plane')
    ax2.scatter(pos_array[:,0]/mks.au,pos_array[:,1]/mks.au,color='black',marker='.',label='Observer Plane')
    ax2.set_xlabel('X-Coordinate (au)',fontsize=15)
    ax2.set_ylabel('Y-Coordinate (au)',fontsize=15)
    ax2.tick_params(axis='both',direction='in')
    ax2.tick_params(axis='both',which='minor',direction='in')
    ax2.tick_params(top=True,right=True)
    ax2.tick_params(which='minor',top=True,right=True)
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax2.legend()

    #This plots the X(t) and Y(t) just to show the cos/sin offset.
    fig3,ax3 = plt.subplots(figsize=(10,8))
    ax3.scatter(t_arr/86400,pos_array[:,0]/mks.au,color='black',marker='x',label='X-coord')
    ax3.scatter(t_arr/86400,pos_array[:,1]/mks.au,color='blue',marker='.',label='Y-coord')
    ax3.set_xlabel('Time (days)',fontsize=15)
    ax3.set_ylabel('Position (au)',fontsize=15)
    ax3.tick_params(axis='both',direction='in')
    ax3.tick_params(axis='both',which='minor',direction='in')
    ax3.tick_params(top=True,right=True)
    ax3.tick_params(which='minor',top=True,right=True)
    ax3.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax3.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax3.legend()

    #This plots the Vz(t) to show there is a non-zero Z-component velocity
    #in the observer reference frame (which would show up as a Vr for the
    #observer on Earth)
    fig4,ax4 = plt.subplots(figsize=(10,8))
    ax4.scatter(t_arr/86400,vel_array[:,2]/1000,color='black',marker='.')
    ax4.set_xlabel('Time (days)',fontsize=15)
    ax4.set_ylabel(r'$\textrm{V}_{\textrm{z}}$ (km s$^{-1}$)',fontsize=15)
    ax4.tick_params(axis='both',direction='in')
    ax4.tick_params(axis='both',which='minor',direction='in')
    ax4.tick_params(top=True,right=True)
    ax4.tick_params(which='minor',top=True,right=True)
    ax4.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax4.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    plt.show()

#This function generates the plots asked for in problem 3. I made up angles
#for i,w,O to iterate through using 3d arrays.
def prob3_plotter(R0,V0):
    #Generate all of the variable values for i,w,O.
    var_params = np.zeros(5,dtype=[('INC','f8'),('wom','f8'),('Om','f8')])
    var_params['INC'] = np.array([20,50,100,150,300])*np.pi/180
    var_params['wom'] = np.array([50,100,150,225,300])*np.pi/180
    var_params['Om'] = np.array([20,100,150,225,300])*np.pi/180

    #Create blank 3d arrays for X,Y position for variable i.
    pos_varINC_array = np.zeros((num_steps,3,5),dtype='f8')
    vel_varINC_array = np.zeros((num_steps,3,5),dtype='f8')
    #Initialize all 3 figures up front.
    fig1,ax1 = plt.subplots(figsize=(10,8))
    fig2,ax2 = plt.subplots(figsize=(10,8))
    fig3,ax3 = plt.subplots(figsize=(10,8))
    #This iterates through all 5 values of i using the same input Parameters
    #from the test_point generation function. e cannot be exactly zero as the
    #root-finding functions are undefined in this case. So I made it really small
    #but still not 0.
    for j in range(5):
        back_half_check = 0
        pos_varINC_array[0,:,j] = R0
        vel_varINC_array[0,:,j] = V0
        a_temp = 1.25*mks.au
        e_temp = 1e-8 #The code requires an non-zero eccentricity. So I used a small value.
        O_temp = 50*np.pi/180
        i_temp = var_params['INC'][j]
        w_temp = 20*np.pi/180
        Ep = 0
        #This part iterates through the X,Y for this loops i-value.
        for k in range(1,num_steps):
            t = t_arr[k]
            if k > (num_steps/2):
                back_half_check = 1
            pos_vec, vel_vec, orb_pos,E_n = pos_vel(a_temp,e_temp,O_temp,i_temp,w_temp,t,T,back_half_check,mag_h,Ep)
            Ep = E_n
            pos_varINC_array[k,:,j],vel_varINC_array[k,:,j] = pos_vec[:],vel_vec[:]
        #Plot this X,Y for the given i. Label them accordinginly.
        plab = 'i = %6.2f$^{\circ}$'%(i_temp*180/np.pi)
        ax1.scatter(pos_varINC_array[:,0,j]/mks.au,pos_varINC_array[:,1,j]/mks.au,marker='.',label=plab)
    ax1.set_xlabel('X-Coordinate (au)',fontsize=15)
    ax1.set_ylabel('Y-Coordinate (au)',fontsize=15)
    ax1.tick_params(axis='both',direction='in')
    ax1.tick_params(axis='both',which='minor',direction='in')
    ax1.tick_params(top=True,right=True)
    ax1.tick_params(which='minor',top=True,right=True)
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax1.legend()

    #This bit will iterate through the values of w.
    pos_varWOM_array = np.zeros((num_steps,3,5),dtype='f8')
    vel_varWOM_array = np.zeros((num_steps,3,5),dtype='f8')
    for j in range(5):
        back_half_check = 0
        pos_varWOM_array[0,:,j] = R0
        vel_varWOM_array[0,:,j] = V0
        a_temp = 1.25*mks.au
        e_temp = 0.5
        O_temp = 50*np.pi/180
        i_temp = 60*np.pi/180
        w_temp = var_params['wom'][j]
        Ep = 0
        #Generate all X,Y for the given w in this iteration.
        for k in range(1,num_steps):
            t = t_arr[k]
            if k > (num_steps/2):
                back_half_check = 1
            pos_vecw, vel_vecw,orb_posw,E_n = pos_vel(a_temp,e_temp,O_temp,i_temp,w_temp,t,T,back_half_check,mag_h,Ep)
            Ep = E_n
            pos_varWOM_array[k,:,j],vel_varWOM_array[k,:,j] = pos_vecw[:],vel_vecw[:]
        #And plot each with the right labels.
        plabWOM = r'$\omega$ = %6.2f$^{\circ}$'%(w_temp*180/np.pi)
        ax2.scatter(pos_varWOM_array[:,0,j]/mks.au,pos_varWOM_array[:,1,j]/mks.au,marker='.',label=plabWOM)
    ax2.set_xlabel('X-Coordinate (au)',fontsize=15)
    ax2.set_ylabel('Y-Coordinate (au)',fontsize=15)
    ax2.tick_params(axis='both',direction='in')
    ax2.tick_params(axis='both',which='minor',direction='in')
    ax2.tick_params(top=True,right=True)
    ax2.tick_params(which='minor',top=True,right=True)
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax2.legend()

    #Finally generate the plot for variable Omega, given a w = pi/2
    pos_varOM_array = np.zeros((num_steps,3,5),dtype='f8')
    vel_varOM_array = np.zeros((num_steps,3,5),dtype='f8')
    for j in range(5):
        back_half_check = 0
        pos_varOM_array[0,:,j] = R0
        vel_varOM_array[0,:,j] = V0
        a_temp = 1.25*mks.au
        e_temp = 0.5
        O_temp = var_params['Om'][j]
        i_temp = 60*np.pi/180
        w_temp = np.pi/2
        Ep = 0
        #This generates all X,Y for a given O.
        for k in range(1,num_steps):
            t = t_arr[k]
            if k > (num_steps/2):
                back_half_check = 1
            pos_vecO, vel_vecO,orb_posO,E_n = pos_vel(a_temp,e_temp,O_temp,i_temp,w_temp,t,T,back_half_check,mag_h,Ep)
            Ep = E_n
            pos_varOM_array[k,:,j],vel_varOM_array[k,:,j] = pos_vecO[:],vel_vecO[:]
        #And plot with the right labels.
        plabOM = r'$\Omega$ = %6.2f$^{\circ}$'%(O_temp*180/np.pi)
        ax3.scatter(pos_varOM_array[:,0,j]/mks.au,pos_varOM_array[:,1,j]/mks.au,marker='.',label=plabOM)
    ax3.set_xlabel('X-Coordinate (m)',fontsize=15)
    ax3.set_ylabel('Y-Coordinate (m)',fontsize=15)
    ax3.tick_params(axis='both',direction='in')
    ax3.tick_params(axis='both',which='minor',direction='in')
    ax3.tick_params(top=True,right=True)
    ax3.tick_params(which='minor',top=True,right=True)
    ax3.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax3.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax3.legend(loc='upper right')

    plt.show()

#Here's the main function for calling all of this. If you want a self-consistent
#test point, select option 1 in the first question. THERE IS NO ERROR-CORRECTION
#IF YOU SELECT OPTION 2 AND GIVE IT GARBAGE. R0,V0 need to be self-consistent
#in option 2 or everything may break and I'm not responsible.
if __name__=='__main__':
    print('Would you like to generate a test point? Or usee your own Pos/Vel Vectors?')
    print('[1] Test Point\n[2] Own Parameters')
    tpoint_select = input('Selection: ')
    #If you want the test_point I selected, here.
    if tpoint_select == '1':
        r0,R0,v0,V0 = gen_test_point(1.25,0.02,60,20,50)
    #This will do everything from a user-input position and velocity, but it better
    #be accurate and consistent. Garbage output means input was garbage, not my
    #fault.
    elif tpoint_select == '2':
        Xus,Yus,Zus = input('Please input X,Y,Z (in cm): ').split(',')
        Vxs,Vys,Vzs = input('Please input Vx,Vy,Vz (in km/s): ').split(',')
        R0 = np.array([np.float(Xus)/100,np.float(Yus)/100,np.float(Zus)/100])
        V0 = np.array([np.float(Vxs)*1000,np.float(Vys)*1000,np.float(Vzs)*1000])
        Rmag0 = vmag(R0)
        r0 = np.array([Rmag0,0,0])
    #This will print the orbital parameters used, whether back-solved from the
    #test point (to make sure the code functions correctly), or generated from
    #the user input.
    print('\nGenerating Orbital Paramters (a,e,i,Omega,w,T,h)')
    print('--------------------------------------------------')
    a,e,O,i,w,T,mag_h = orbit_params(R0,V0)
    #Now we iterate through time to build the array for problem 2.
    num_steps = 1000
    step_size = T / num_steps
    #Initialize full XYZ, VxVyVz, and xyz arrays.
    pos_array = np.zeros((num_steps,3),dtype='f8')
    vel_array = np.zeros((num_steps,3),dtype='f8')
    orb_array = np.zeros((num_steps,3),dtype='f8')
    #Set the starting points.
    pos_array[0,:] = R0
    vel_array[0,:] = V0
    orb_array[0,:] = r0
    back_half_check = 0
    #time as an array so I know that XYZ steps use the same time values that the
    #plots will use.
    t_arr = np.zeros(num_steps,dtype='f8')
    t_arr = np.array([j*step_size for j in range(num_steps)])
    Ep = 0 #First angle guess, since it starts at 0.
    for j in range(1,num_steps):
        t = t_arr[j]
        if j > (num_steps/2):
            back_half_check = 1
        pos_vec, vel_vec,orb_pos,E_n = pos_vel(a,e,O,i,w,t,T,back_half_check,mag_h,Ep)
        pos_array[j,:],vel_array[j,:],orb_array[j,:] = pos_vec,vel_vec,orb_pos
        Ep = E_n #this allows the next iteration to have a better starting point
        #The program will not converge for some angle is 0 is always used.
    #Ask the user if they want to plot stuff. Problem 2 generates 3 plots.
    #Problem 3 generates 15 plots on 3 figures.
    p2_choice = input('\nWould you like to generate the plots for problem 2? [y/n]: ')
    if ((p2_choice == 'y')|(p2_choice=='Y')|(p2_choice=='yes')|(p2_choice=='YES')):
        prob2_plotter(pos_array,orb_array,vel_array,t_arr)
    p3_choice = input('\nWould you like to generate the plots for problem 3? [y/n]: ')
    if ((p3_choice == 'y')|(p3_choice=='Y')|(p3_choice=='yes')|(p3_choice=='YES')):
        prob3_plotter(R0,V0)
    #And we're done. No active menu this time.
    print('\nGoodbye!')
