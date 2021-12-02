
import matplotlib.pylab as plt
import numpy as np
import scipy.integrate as sp
import scipy.special as sp1

'---------------------------'
'--COSMOLOGICAL PARAMETERS--'
'-------(Planck 2015)-------'
'---------------------------'

Omega_b_h2 = 0.0223         # Physical baryon density parameter
Omega_c_h2 = 0.1188         # Physical dark matter density parameter
W  = -1                     # Equation state of dark energy 
Omega_b = 0.0486            # Baryon density parameter
Omega_c = 0.2589            # Dark Matter density parameter
Omega_m = 0.3089            # Matter density parameter
Omega_lambda = 0.6911       # Dark Energy density parameter
h = 0.6774
H0 = (100*h)*(1e3/3.086e22) # Hubble constant in SI units

"NOTE:- Omega_k = 0 for LCDM model"
"NOTE:- Omega_rad = 0 if we consider the above values"

'---------------------------'
'------OTHER CONSTANTS------'
'---------------------------'

hp = 6.626e-34              # Planck constant
c = 3.0e8                   # Speed of light
me = 9.1e-31                # Mass of electron
kB = 1.38e-23               # Boltzmann constant
sigma_T = 6.6524e-29        # Thomson scattering cross-section 
sigma_SB = 5.67e-8          # Stefan-Boltzmann constant
T0 = 2.725                  # CMB Temperature at z = 0
del_E = 1.632e-18           # Ly_alpha transition energy
E2 = 5.44e-19               # Energy of 2s orbital of hydrogen
E1 = 2.176e-18              # Energy of 1s orbital of hydrogen
G = 6.67e-11                # Gravitational constant

'-------------------------------------------------------------'
'-------CONSTANTS FOR COLLISIONAL COUPLING COEFFEICIENTS------'
'-------------------------------------------------------------'

T_star = 0.068             
A_10 = 2.85e-15          # Einstein coefficient for spontaneous emission

'-----------------------------------------'
'''FUNCTION TO DETERMINE CMB TEMPERATURE'''
'-----------------------------------------'

def T_CMB(red):
    T_gamma = T0*(1.0 + red)
    return T_gamma

'------------------------------------------'
'''FUNCTION TO DETERMINE HUBBLE PARAMETER'''
'------------------------------------------'

def H(red):
    H_z= H0*(Omega_m*(1.0 + red)**3.0)**0.5
    return H_z

'---------------------------------------------------------'
'''FUNCTION TO DETERMINE NEUTRAL HYDROGEN NUMBER DENSITY'''
'---------------------------------------------------------'

def nH(red):
    NH = 8.403*Omega_b_h2*(1.0 + red)**3.0
    return NH

'-----------------------------------------------------------'
'''FUNCTION TO DETERMINE RECOMBINATION COEFFICIENT (alpha)'''
'-----------------------------------------------------------'

def alpha_e(Tg):
    a = 4.309
    b = -0.6166
    cp = 0.6703
    d = 0.5300
    F = 1.14             # Fudge factor
    t = Tg/(1.0e4)
    
    alpha = F*1.0e-19*((a*t**b)/(1.0 + cp*t**d))
    return alpha

'------------------------------------------------------------'
'''FUNCTION TO DETERMINE PHOTOIONIZATION COEFFICIENT (beta)'''  
'------------------------------------------------------------'
 
def beta_e(Tg):
    beta = alpha_e(Tg)*2.4093e21*Tg**(1.5)*np.exp(-39420.289/Tg)
    return beta

'---------------------------'
'''FUNCTION TO DETERMINE C'''
'---------------------------'

def C1(red,x, Tg):
    K = 7.1898e-23/(H(red))
    Lambda = 8.22458 
    
    Cr = (1.0 + K*Lambda*(1.0 - x)*nH(red))/(1.0 + K*(Lambda + beta_e(Tg))*(1.0 - x)*nH(red))
    return Cr

'--------------------------------------'
'''FUNCTION TO DETERMINE BOOST FACTOR'''
'--------------------------------------'

def B(red):
    B_z = 1.0 + (((1.6e5)*sp1.erfc((1.0 + red)/20.5))/((1.0 + red)**(1.54)))
    return B_z

'---------------------------------------------'
'''FUNCTION TO DETERMINE DARK MATTER DENSITY'''
'---------------------------------------------'

def rho_DM(red):
    rho_DM_eng =  (Omega_m-Omega_b)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_DM_eng

'----------------------------------------'
'''FUNCTION TO DETERMINE MATTER DENSITY'''
'----------------------------------------'

def rho_M(red):
    rho_M_eng = (Omega_m)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_M_eng    

'----------------------------------------'
'''FUNCTION TO DETERMINE BARYON DENSITY'''
'----------------------------------------'

def rho_B(red):
    rho_B_eng = (Omega_b)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_B_eng

'------------------------------'
'''FUNCTION TO DETERMINE u_th'''
'------------------------------'

def u_th(Tg, Tx, dm):
    mx = dm
    return c*2.936835e-7*(((Tg/mb) + (Tx/mx))**(0.5))

'------------------------------'
'''FUNCTION TO DETERMINE F(r)'''
'------------------------------'
  
def F_r(vel_xb, Tg, Tx, dm):
    u_therm = u_th(Tg, Tx, dm)
    rv = vel_xb/u_therm
    F = sp1.erf(rv/np.sqrt(2.0)) - np.sqrt(2.0/np.pi)*np.exp((-rv**2.0)/2.0)*rv
    return F

'----------------------------------------'
'''FUNCTION TO DETERMINE F(r)/Vel_xb^2'''
'----------------------------------------'


def Fr_by_velxb2(vel_xb, Tg, Tx, dm):
    u_therm = u_th(Tg, Tx, dm)
    rv = vel_xb/u_therm
    if rv >= 0.1:
        F = sp1.erf(rv/np.sqrt(2.0)) - np.sqrt(2.0/np.pi)*np.exp((-rv**2.0)/2.0)*rv
        F = F/vel_xb**2
        return F
    else:
        F = np.sqrt(2.0/np.pi)*(rv/3.0 - rv**3.0/10.0 + rv**5.0/56.0)
        F = F/u_therm**2
        return F
        
'----------------------------------------'
'''FUNCTION TO DETERMINE F(r)/Vel_xb'''
'----------------------------------------'


def Fr_by_velxb(vel_xb, Tg, Tx, dm):
    u_therm = u_th(Tg, Tx, dm)
    rv = vel_xb/u_therm
    if rv >= 0.1:
        F = sp1.erf(rv/np.sqrt(2.0)) - np.sqrt(2.0/np.pi)*np.exp((-rv**2.0)/2.0)*rv
        F = F/vel_xb
        return F
    else:
        F = np.sqrt(2.0/np.pi)*(rv**2.0/3.0 - rv**4.0/10.0 + rv**6.0/56.0)
        F = F/u_therm
        return F


'---------------------------------------'
'''FUNCTION TO DETERMINE THE DRAG TERM'''
'---------------------------------------'

def Drag_vxb(vel_xb, red, Tg, Tx, sigma_45, dm):
    mx = dm
    D_Vxb2 = 1.*2.63e+7*h*(Omega_m**0.5)*((1.0 + red)**0.5)*sigma_45*Fr_by_velxb2(vel_xb, Tg, Tx, dm)/(mb + mx)
    return D_Vxb2

'---------------------------------------'
'''FUNCTION TO DETERMINE Q_b_coupling'''
'---------------------------------------'

def Q_b_coupling(Tx, Tg, red, vel_xb, sigma_45, dm):
    mx = dm
    u_therm = u_th(Tg, Tx, dm)
    rv = vel_xb/u_therm
    Q_b1 = (2.0/3.0)*2.10e+7*(Omega_m - Omega_b)*(h**2.0)*((1.0 + red)**0.5)*sigma_45*np.exp(-(rv**2)/2.0)*(mb/((mb + mx)**2.0))*(Tx - Tg)/(h*(Omega_m**0.5)*u_therm**3)/1.0
    return Q_b1


'-----------------------------'
'''FUNCTION TO DETERMINE Q_b_drag'''
'-----------------------------'

def Q_b_drag(Tx, Tg, red, vel_xb, sigma_45, dm):
    mx = dm
    Q_b_d = 1.*2.26e+3*(Omega_m - Omega_b)*h**2*(1.0+red)**0.5*sigma_45*(mb*mx/((mb + mx)**2.0))*Fr_by_velxb(vel_xb, Tg, Tx, dm)/(h*Omega_m**0.5) 
    return Q_b_d

'-----------------------------'
'''FUNCTION TO DETERMINE Q_x_coupling'''
'-----------------------------'

def Q_x_coupling(Tx, Tg, red, vel_xb, sigma_45, dm):
    mx = dm
    u_therm = u_th(Tg, Tx, dm)
    rv = vel_xb/u_therm
    Q_x1 = (2.0/3.0)*2.10e+7*Omega_b*(h**2.0)*((1.0 + red)**0.5)*sigma_45*np.exp(-(rv**2)/2.0)*(mx/((mb + mx)**2.0))*(Tg - Tx)/(h*(Omega_m**0.5)*u_therm**3)/1.0
    return Q_x1


'-----------------------------'
'''FUNCTION TO DETERMINE Q_x_drag'''
'-----------------------------'

def Q_x_drag(Tx, Tg, red, vel_xb, sigma_45, dm):
    mx = dm
    Q_x_d = 1.*2.26e+3*Omega_b*h**2*(1.0+red)**0.5*sigma_45*(mb*mx/((mb + mx)**2.0))*Fr_by_velxb(vel_xb, Tg, Tx, dm)/(h*Omega_m**0.5)  
    return Q_x_d

    

'----------------------------------------------------------------------------------'
'''FUNCTION TO DETERMINE GAS KINETIC TEMPERATURE (Tg) AND IONIZATION FRACTION (x)'''
'----------------------------------------------------------------------------------'
    
def func(r,red, cross_in, cross_an, dm):
    Tg = r[0]
    x = r[1]
    Tx = r[2]
    V_xb = r[3]
    
    
    f_Tg = ((2.0*Tg)/(1.0 + red)) - ((2.70877e-20*(T_CMB(red) - Tg)*(1.0 + red)**(1.5)*(x/(1.0 + x)))/(H0*np.sqrt(Omega_m))) - Q_b_coupling(Tx, Tg, red, V_xb, cross_in, dm)-Q_b_drag(Tx, Tg, red, V_xb, cross_in, dm)/1.0 - 1.056e31*(((Omega_m - Omega_b)**2)/(Omega_b*(Omega_m**0.5)))*h*((1.0 + red)**0.5)*B(red)*((1.0 + 2*x)/(1.0 + x))*(cross_an/dm)
    f_x =  (C1(red,x,Tg)*(alpha_e(Tg)*x**2*nH(red) - beta_e(T_CMB(red))*(1.0 - x)*np.exp(-118260.87/Tg)))/(H(red)*(1.0 + red)) - 2.1875e8*(((Omega_m - Omega_b)**2)/(Omega_b*(Omega_m**0.5)))*h*((1.0 + red)**0.5)*((1.0/E1) + ((1.0 - C1(red,x,Tg))/del_E))*(1.0 - x)*B(red)*(cross_an/dm)
    f_Tx = ((2.0*Tx)/(1.0 + red)) - Q_x_coupling(Tx, Tg, red, V_xb, cross_in, dm) -Q_x_drag(Tx, Tg, red, V_xb, cross_in, dm)/1.0
    f_Vxb = (V_xb/(1.0 + red)) + Drag_vxb(V_xb, red, Tg, Tx, cross_in, dm)
    
    return np.array([f_Tg, f_x, f_Tx, f_Vxb], float)



'------------------------'
'''WITHOUT ANNIHILATION'''
'------------------------'

def func2(r2,red, cross_in, dm):
    
    Tg2 = r2[0]
    x2 = r2[1]
    Tx2 = r2[2]
    V_xb2 = r2[3]
    
    
    f_Tg2 = ((2.0*Tg2)/(1.0 + red)) - ((2.70877e-20*(T_CMB(red) - Tg2)*(1.0 + red)**(1.5)*(x2/(1.0 + x2)))/(H0*np.sqrt(Omega_m))) - Q_b_coupling(Tx2, Tg2, red, V_xb2, cross_in, dm)-Q_b_drag(Tx2, Tg2, red, V_xb2, cross_in, dm)/1.0
    f_x2 =  (C1(red,x2,Tg2)*(alpha_e(Tg2)*x2**2*nH(red) - beta_e(T_CMB(red))*(1.0 - x2)*np.exp(-118260.87/Tg2)))/(H(red)*(1.0 + red)) 
    f_Tx2 = ((2.0*Tx2)/(1.0 + red)) - Q_x_coupling(Tx2, Tg2, red, V_xb2, cross_in, dm) -Q_x_drag(Tx2, Tg2, red, V_xb2, cross_in, dm)/1.0
    f_Vxb2 = (V_xb2/(1.0 + red)) + Drag_vxb(V_xb2, red, Tg2, Tx2, cross_in, dm)
    
    return np.array([f_Tg2, f_x2, f_Tx2, f_Vxb2], float)

def prob_func(v_xb):
    
    return 4*np.pi*v_xb**2*(np.exp((-3*v_xb**2)/(2*Vrms**2)))*(1.0/((2.0/3)*np.pi*Vrms**2)**(1.5))
    

zf = 9.0                    # Final redshift
zi = 1010.0                  # Initial redshift
del_z = -0.01       # Step-size

z = np.arange(zi, zf, del_z)


Vxb = np.arange(2.0e-6*c, 5.0e-4*c, 2.0e-5*c)
Vrms = 1e-4*c



'----------CMB TEMPERATURE---------'

T_gamma = T0*(1.0 + z)

T_gas_xb = []
T_gas2_xb = []
T_spin_xb = []
T_spin2_xb = []
T_21_xb = []
T_21_2_xb = []
T_dark_xb = []
T_dark2_xb = []
x_points_xb = []
x_points2_xb = []
P_Vxb = []

for vxb in Vxb:
    
    '----------------------------------------------------'
    '''---------------INITIAL CONDITIONS---------------'''
    '---- (Tg = T_CMB at z = 1010 and x = 0.05497369) ----'
    '----------------------------------------------------'

    r0 = np.array([T_CMB(zi), 0.05497369, 0.0, vxb], float)    # Initial conditions for Tg and x

    
    
    '''------SOLVING THE EQUATIONS--------'''

    mx = 0.0001                          # Mass of dark matter particle
    mb = 0.938                        # Mass of baryonic particle 
    sigma_45 = 1000.0                   # Interaction cross-section
    sigma_v =  (1e-37)            # Annihilation cross-section    
    
    
    r = sp.odeint(func, r0, z, args = (sigma_45, sigma_v, mx))            # Solving the coupled differential equation
    T_gas = r[:,0]                        # Stores Tg as an array
    x_points = r[:,1]                     # Stores x as an array
    T_dark = r[:,2]                       # Stores Tx as an array
    Vel_xb = r[:,3]                       # Stores V_xb as an array

    r2 = sp.odeint(func2, r0, z, args = (sigma_45, mx))
    T_gas2 = r2[:,0]                      # Stores Tg as an array
    x_points2 = r2[:,1]                     # Stores x as an array
    T_dark2 = r2[:,2]                       # Stores Tx as an array
    Vel_xb2 = r2[:,3]                       # Stores V_xb as an array
    
    
    

    '-------CALCULATION OF COLLISIONAL COUPLING COEFFICIENT------'

    K_HH = 3.1e-17*T_gas**(0.357)*np.exp(-32.0/T_gas)     # Using the fitting formula
    K_HH2 = 3.1e-17*T_gas2**(0.357)*np.exp(-32.0/T_gas2)
    nHI = 8.403*Omega_b_h2*(1.0 + z)**3.0

    C_10 = K_HH*nHI
    C_10_2 = K_HH2*nHI

    x_c = (T_star*C_10)/(A_10*T_gamma)                   # Collisional coupling coefficient 
    x_c2 = (T_star*C_10_2)/(A_10*T_gamma)
    
    
    

    '------CALCULATION OF SPIN TEMPERATURE------'

    T_spin = ((1.0 +  x_c)*T_gas*T_gamma)/(x_c*T_gamma + T_gas)
    T_spin2 = ((1.0 +  x_c2)*T_gas2*T_gamma)/(x_c2*T_gamma + T_gas2)
    
    T_21 = 23*((0.15/Omega_m)*((1.0 + z)/10))**(0.5)*(Omega_b*h/0.02)*(1.0 - (T_gamma/T_gas))
    T_21_2 = 23*((0.15/Omega_m)*((1.0 + z)/10))**(0.5)*(Omega_b*h/0.02)*(1.0 - (T_gamma/T_gas2))
    
    
    
    
    P_Vxb.append(prob_func(vxb))
    T_gas_xb.append(T_gas)
    T_gas2_xb.append(T_gas2)
    T_spin_xb.append(T_spin)
    T_spin2_xb.append(T_spin2)
    T_21_xb.append(T_21)
    T_21_2_xb.append(T_21_2)
    T_dark_xb.append(T_dark)
    T_dark2_xb.append(T_dark2)
    x_points_xb.append(x_points)
    x_points2_xb.append(x_points2)
    
    
    
    
T_21_avg = []
T_gas_avg = []
T_dark_avg = []
x_points_avg = []
T_spin_avg = []

T_21_2_avg = []
T_gas2_avg = []
T_dark2_avg = []
x_points2_avg = []
T_spin2_avg = []


for o in range(len(P_Vxb)):
    T_21_avg.append(T_21_xb[o]*P_Vxb[o])
    x_points_avg.append(x_points_xb[o]*P_Vxb[o])
    T_gas_avg.append(T_gas_xb[o]*P_Vxb[o])
    T_dark_avg.append(T_dark_xb[o]*P_Vxb[o])
    T_spin_avg.append(T_spin_xb[o]*P_Vxb[o])
    
    T_21_2_avg.append(T_21_2_xb[o]*P_Vxb[o])
    x_points2_avg.append(x_points2_xb[o]*P_Vxb[o])
    T_gas2_avg.append(T_gas2_xb[o]*P_Vxb[o])
    T_dark2_avg.append(T_dark2_xb[o]*P_Vxb[o])
    T_spin2_avg.append(T_spin2_xb[o]*P_Vxb[o])

T_b_avg = sum(T_21_avg)/(1.0*sum(P_Vxb))
x_avg = sum(x_points_avg)/(1.0*sum(P_Vxb))
T_g_avg = sum(T_gas_avg)/(1.0*sum(P_Vxb))
T_s_avg = sum(T_spin_avg)/(1.0*sum(P_Vxb))
T_x_avg = sum(T_dark_avg)/(1.0*sum(P_Vxb))

T_b2_avg = sum(T_21_2_avg)/(1.0*sum(P_Vxb))
x2_avg = sum(x_points2_avg)/(1.0*sum(P_Vxb))
T_g2_avg = sum(T_gas2_avg)/(1.0*sum(P_Vxb))
T_s2_avg = sum(T_spin2_avg)/(1.0*sum(P_Vxb))
T_x2_avg = sum(T_dark2_avg)/(1.0*sum(P_Vxb))


print(T_b_avg[99280], T_b2_avg[99280])

'''Temp_data_with_DM = np.array([z, T_gas, T_spin, T_gamma, T_dark])
Temp_data_with_DM = Temp_data_with_DM.T

Ionize_frac_data_with_DM = np.array([z, x_points])
Ionize_frac_data_with_DM = Ionize_frac_data_with_DM.T'''

#np.savetxt('Various_Temp_vs_z_with_DM_10p0GeV_0.txt', Temp_data_with_DM)
#np.savetxt('Ionization_Frac_vs_z_with_DM_10p0GeV_0.txt', Ionize_frac_data_with_DM)

#plt.loglog(1.0 + z, T_spin, 'r', label = 'T_spin_AN')
#plt.loglog(1.0 + z, T_spin2, 'b', label = 'T_spin_no_AN')
#plt.loglog(1.0 + z, T_gamma, 'g--', label = 'T_cmb')
'''plt.loglog(1.0 + z, T_g_avg, 'r', label = r'$f^2_{DM}f_{eff}<\sigma v>$ = ')
plt.loglog(1.0 + z, T_g2_avg, 'b', label = r'$f^2_{DM}f_{eff}<\sigma v>$ = 0')
plt.loglog(1.0 + z, T_x_avg, 'r--', label = r'$f^2_{DM}f_{eff}<\sigma v>$ = ')
plt.loglog(1.0 + z, T_x2_avg, 'b--', label = r'$f^2_{DM}f_{eff}<\sigma v>$ = 0')
#plt.loglog(1.0 + z, T_dark)
plt.xlim(10.0, 1000.0)
#plt.ylim(1e-2, 1e3)
plt.legend(loc = 'best')'''

plt.plot(1.0 + z, T_b_avg, 'r', label = r'$f^2_{DM}f_{eff}<\sigma v>$ = ')
plt.plot(1.0 + z, T_b2_avg, 'b', label = r'$f^2_{DM}f_{eff}<\sigma v>$ = 0')
plt.xlim(10.0, 300.0)
#plt.ylim(-60.0, 0.0)
plt.legend(loc = 'best')

'''c_s = c*np.sqrt(3*kB*T_gas/(mb*1.6e-10))
plt.loglog(1.0 + z, (Vel_xb/(1e3*(1.0 + z))))
plt.loglog(1.0 + z, (c_s/(1e3*(1.0 + z))))
plt.xlim(30.0, 250.0)
plt.ylim(0.01, 0.035)'''

plt.show()
 





