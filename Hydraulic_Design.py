###Imports###
import numpy as np
import math


###Variables for hydraulic design###
N = 13 #Number of tubes
L = 0.35 #Length
Y = 14e-3 #Tube pitch
D_sh = 64e-3 #Acrylic shell diameter
mu = 6.51e-4 #Dynamic viscosity (kg/ms)
do = 8e-3 #Tube outer diameter
di = 6e-3 #Tube inner diameter
N_baffle = 9 #Number of baffles
d_noz = 19e-3 #THIS NEEDS CHANGING AFTER CONSULTATION ABOUT 19MM USED INSTEAD OF 24.5????????
rho = 990.1 #kg/m3 (water density at 40C)
a = 0.34 #Use 0.2 for triangular pitch. Use 0.34 for square pitch. Constant used for shell pressure drop, see handout top of page 4
m_dot_h = 0.45 #kg/s (total mass flow of hot stream (not per tube))
m_dot_c = 0.50 #kg/s
kc = 0.45 #For testing
ke = 0.8 #For testing



###Definitions of correlations###
B = L/(N_baffle+1) #Baffle spacing
A_sh = D_sh*(Y-do)*B/Y #This correlation is only approx (see notes eqn 6)
V_sh = m_dot_c/(rho*A_sh) #Measure of shell velocity


###Functions for calculations###
def Dp_in_out(v_tube,kc,ke):
    #Returns entrance and exit pressure losses given kc and ke. These are found from edge_coefficients
    return 0.5*rho*v_tube**2*(kc+ke)    
def v_tube(m_dot_h, N):
    #Returns the tube veloctiy for the given input parameters
    if (m_dot_h/N) /(np.pi*rho*0.25*di**2) < 0.393:
        print("Caution, equation (4) in datasheet not valid as tube Reynolds too low")
    else:
        return (m_dot_h/N) /(np.pi*rho*0.25*di**2)
def sigma(N):
    #Returns sigma used for ke and kc evalutation
    return (N*0.25*np.pi*di**2)/(0.25*np.pi*D_sh**2) #Free area to total area
def Dp_tube(f,v_tube):
    #Returns pressure loss due to friction in tube
    return (f*L*rho*v_tube**2)/(di*2)
def friction_factor(v_tube):
    #From worked example top of page 2. we can find friction factor. Returns friction factor, f
    return (1.82*math.log10(rho*v_tube*di/mu)-1.64)**(-2)
def v_noz(m_dot,dia):
    #Returns the flow speed for a given mass flow rate and diameter. Use v_tube for hot side velocity to account for N tubes
    return m_dot/(rho*0.25*np.pi*dia**2)
def Dp_shell(V_sh,N):
    #Returns the pressure loss from the dubious empirical relation given on handout page 3 at bottom
    return 4*a*(V_sh*do*rho/mu)**-0.15*N*rho*V_sh**2


def edge_coefficients(sigma,):

    return kc,ke



v_tube = v_tube(m_dot_h,N)
f = friction_factor(v_tube)
v_noz_h = v_noz(m_dot_h,d_noz)
v_noz_c = v_noz(m_dot_c,d_noz)

Dp_h_overall = Dp_tube(f,v_tube) + Dp_in_out(v_tube, kc, ke) + rho*v_noz_h**2 #Last term accounts for BOTH nozzles
Dp_c_overall = Dp_shell(V_sh,N) + rho*v_noz_c**2 #Last term acoounts for BOTH nozzles
print(Dp_shell(V_sh,N))