###Imports###
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import interpolate
from intersect import intersection #Used for finding mass flow rates by iteration


###Figure 6 Dataset###
dataset_Qc = 1e-3*np.array([0.5833,0.5083,0.4750,0.4250,0.3792,0.3417,0.2958,0.2583,0.2125,0.1708]) #m3/s
dataset_Dpc = 1e5*np.array([0.1113,0.2157,0.2538,0.3168,0.3613,0.4031,0.4511,0.4846,0.5181,0.5573]) #in Pa not bar
dataset_Qh = 1e-3*np.array([0.4583,0.4236,0.4010,0.3611,0.3125,0.2639,0.2222,0.1597,0.1181,0.0694]) #m3/s
dataset_Dph = 1e5*np.array([0.1333,0.1756,0.2024,0.2577,0.3171,0.3633,0.4233,0.4784,0.5330,0.5715]) #in Pa not bar



###Variables for hydraulic design###
N = 14 #Number of tubes total (not per pass)
N_pass_h = 2 #Number of passes of hot flow if multi-pass is used
N_pass_c = 2 #Number of passes of cold flow if multi-pass is used
L = 0.23 #Length of a single tube
Y = 15e-3 #Tube pitch
D_sh = 64e-3 #Acrylic shell diameter
mu = 6.51e-4 #Dynamic viscosity (kg/ms)
do = 8e-3 #Tube outer diameter
di = 6e-3 #Tube inner diameter
N_baffle = 12 #Number of baffles PER SIDE, CAREFUL. Each baffle often on both sides! (multipass shell)
d_noz = 25e-3 #THIS NEEDS CHANGING AFTER CONSULTATION ABOUT 20MM USED INSTEAD OF 24.5????????
rho = 1/0.001008 #kg/m3 (water density at 40C)
a = 0.34 #Use 0.20 for triangular pitch. Use 0.34 for square pitch. Constant used for shell pressure drop, see handout top of page 4
L_plenum = 41e-3 #Plenum length


###Definitions of correlations###
B = L/(N_baffle+1) #Baffle spacing PER SIDE (multipass shell)

#------------------
A_sh = D_sh*(Y-do)*B/(Y*N_pass_c) #This correlation is only approx (see notes eqn 6)         #why divide by N_pass_c????
sigma = N*0.25*np.pi*di**2/(0.25*np.pi*D_sh**2) #Free area to total area. Is this equation right for multipass? i.e sigma does not vary with multipass vs signle pass?
sigma_plenum = (np.pi*0.25*d_noz**2)/(0.5*(d_noz+D_sh)*L_plenum) #Area ratio for nozzle into plenum. Assume plenum is a rectangle area L by (d_noz+D_sh)*0.5            #why need this?????
#----------------------


###Functions for calculations###
def Dp_in_out(v_tube,kc,ke):
    #Returns entrance and exit pressure losses given kc and ke. These are found from edge_coefficients. Is this equation right for multipass?
    return 0.5*rho*v_tube**2*(kc+ke)*N_pass_h
def v_tube(m_dot, N):
    #Returns the tube veloctiy for the given input parameters taking into account multi-pass
    if (m_dot*N_pass_h/N)/(np.pi*rho*0.25*di**2) < 0.393:
        print("Caution, equation (4) in datasheet not valid as tube Reynolds too low")
    else:
        return (m_dot*N_pass_h/N)/(np.pi*rho*0.25*di**2)
def Dp_tubes(f,v_tube):
    #Returns pressure loss due to friction in tubes for all passes (if multi-pass)
    return (f*L*N_pass_h*rho*v_tube**2)/(di*2)
def friction_factor(v_tube):
    #From worked example top of page 2. we can find friction factor. Returns friction factor, f
    return (1.82*math.log10(rho*v_tube*di/mu)-1.64)**(-2)
def v_noz(m_dot,dia):
    #Returns the flow speed for a given mass flow rate and diameter. Use v_tube for hot side velocity to account for N tubes
    return m_dot/(rho*0.25*np.pi*dia**2)
def Dp_shell(V_sh,N):
    #Returns the pressure loss from the dubious empirical relation given on handout page 3 at bottom
    return 4*a*(V_sh*do*rho/mu)**-0.15*(N/N_pass_c)*rho*V_sh**2
def edge_coefficients(sigma):
    #Given a value of sigma, returns kc,ke from Fig 7 (assuming Re=10,000)
    kc = 0.4*(1-sigma) + 0.09 #From Fig 7. using Reynold's as approx 10,000
    ke = sigma**2 - 2.1 * sigma + 1 #From curve fit of Fig 7, again using Re=10,000
    return kc,ke
def plenum_loss(v_noz,kc,ke):
    #Large loss due to mixing in plenum from nozzle entry and viscous mixing. Returns pressure loss. Scaled using Longley 2022 data
    return 9*rho*v_noz**2*(kc+ke)
def m_dot_H(Dp_h_overall):
    #Returns hot mass flow rate for a given pressure drop along hot side. This uses Fig 6 data from handout
    a,b,c,d = np.polyfit(dataset_Dph, dataset_Qh, 3)
    return rho*(a*Dp_h_overall**3+b*Dp_h_overall**2+c*Dp_h_overall+d)
def m_dot_C(Dp_c_overall):
    #Returns cold mass flow rate for a given pressure drop along cold side. This uses Fig 6 data from handout
    a,b,c,d = np.polyfit(dataset_Dpc, dataset_Qc, 3)
    return rho*(a*Dp_c_overall**3+b*Dp_c_overall**2+c*Dp_c_overall+d)



###Iteration to find correct mass flow rates###
m_dot_h_actual = []#Returns flow rate due to guessed flow rate.
m_dot_h_guess = [] #Guess actual flow rate
m_dot_c_actual = []
m_dot_c_guess = []

for m_dot in np.linspace(0.15,0.7,100):
    v_t = v_tube(m_dot,N)
    f = friction_factor(v_t)
    v_noz_h = v_noz(m_dot,d_noz)
    kc1,ke1 = edge_coefficients(sigma)[0],edge_coefficients(sigma)[1] #For tube loss
    kc2,ke2 = edge_coefficients(sigma_plenum)[0],edge_coefficients(sigma_plenum)[1] #For plenum loss
    Dp_h_overall = Dp_tubes(f,v_t) + Dp_in_out(v_t, kc1, ke1) + plenum_loss(v_noz_h,kc2,ke2) + rho*v_noz_h**2 #Last term accounts for BOTH nozzles
    m_dot_h_guess.append(m_dot)
    m_dot_h_actual.append(m_dot_H(Dp_h_overall))

m_dot_h = intersection(m_dot_h_guess,m_dot_h_actual,m_dot_h_guess,m_dot_h_guess)[0][0] #Finds intersection of data i.e actual value of m_dot
print("This configuration has hot mass flow:", round(m_dot_h,5),"kg/s (",round(m_dot_h/rho*1000,3),"ltr/s )")

for m_dot in np.linspace(0.15,0.7,100):
    v_noz_c = v_noz(m_dot,d_noz)
    V_sh = m_dot/(rho*A_sh) #Measure of shell velocity

    #----------------------
    Dp_c_overall = Dp_shell(V_sh,N) + 2.3*N_baffle*N_pass_c*rho*V_sh**2 + rho*v_noz_c**2 #Last term accounts for BOTH nozzles. 2nd term is calibrated by Longley 2022 and 2017B data on Moodle
    #-------------------------------
    
    
    m_dot_c_guess.append(m_dot)
    m_dot_c_actual.append(m_dot_C(Dp_c_overall))
#plt.plot(m_dot_c_guess,m_dot_c_actual)
#plt.plot(m_dot_c_guess,m_dot_c_guess)
#plt.xlabel('Guess')
#plt.ylabel('Actual')
#plt.show()
m_dot_c = intersection(m_dot_c_guess,m_dot_c_actual,m_dot_c_guess,m_dot_c_guess)[0][0] #Finds intersection of data i.e actual value of m_dot
print("This configuration has cold mass flow:", round(m_dot_c,5),"kg/s (",round(m_dot_c/rho*1000,3),"ltr/s )")