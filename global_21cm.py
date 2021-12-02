
#This code just solves the coupled equations for a given mass and cross section. 
#It generates the plots of Tg and Tx, T_spin, and T_21 
#(with averaging over velocity) as given in Munoz's (2015) paper. 

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
Omega_m = 0.3               # Matter density parameter
Omega_lambda = 0.6911       # Dark Energy density parameter
h = 0.6774
H0 = (100*h)*(1e3/3.086e22) # Hubble constant in SI units

"NOTE:- Omega_k = 0 for LCDM model"
"NOTE:- Omega_rad = 0 if we consider the above values"







'--------------------------------------------'
'------OTHER CONSTANTS (all in SI unit)------'
'--------------------------------------------'

hp = 6.626e-34                      # Planck constant
c = 3.0e8                           # Speed of light
me = 9.1e-31                        # Mass of electron
kB = 1.38e-23                       # Boltzmann constant
sigma_T = 6.6524e-29                # Thomson scattering cross-section 
sigma_SB = 5.67e-8                  # Stefan-Boltzmann constant
T0 = 2.725                          # CMB Temperature at z = 0
del_E = 1.632e-18                   # Ly_alpha transition energy
E2 = 5.44e-19                       # Energy of 2s orbital of hydrogen
E1 = 2.176e-18                      # Energy of 1s orbital of hydrogen
G = 6.67e-11                        # Gravitational constant

G_CGS = 6.67e-8                     # Gravitational constant in cm3 gm-1 sec-2
Mpc = 3.085677581e24                # Mega-Parsec in centimetres
H0_CGS = (100*h)*(1e3/3.086e22)     # Hubble constant ! 100*h*10**5.*Mpc**-1., same in SI and CGS unit
kB_CGS = 1.38e-16                   # Boltzmann constant ! (cgs) ! cm^2 g s^-2 K^-1 or ergs/K
c_CGS = 3.0e10                      # Speed of light ! (cgs) 
nu_alpha = 2.47e15                  # Lyman-alpha frequency in Hz
hp_CGS = 6.626e-27                  # Planck constant ! (cgs)





'----------------------------------------------------------'

E1s = 2.176e-11                     # in erg

Xray_ON = float(input('Enter 1 to include Xray photons, 0 to exlude them: ')) # write '1' to enable X-ray photons, '0' to disable it
print('--------------------------------------------------------------------')
Ly_alpha_ON = float(input('Enter 1 to include Ly-alpha photons, 0 to exlude them: ')) # write '1' to enable Lyman-alpha photons, '0' to disable it
print('--------------------------------------------------------------------')



if (Xray_ON == 1 or Ly_alpha_ON == 1):
    
    source_model = float(input('Source model (1 = STARBURST | 2 = SNR | 3 = MINIQUASARS |): '))      # 1 = STARBURST ! 2 = SNR ! 3 = MINIQUASARS !
    
    f_star = float(input('Enter Star formation efficieny: ' ))     
    
    f_star = f_star/0.01                #in units of 0.01    # The data of Xray heating by Raghu is for f_star = 0.01. Change the value accordingly
    
else:

    source_model = 1       # doesn't matter
    f_star = 1.0           # doesn't matter
    
print('--------------------------------------------------------------------')


DM_B_ON = float(input('Enter 1 to include DM-b interaction, 0 to exlude it: ')) # write '1' to enable DM-b interaction, '0' to disable it
print('--------------------------------------------------------------------')


'---------------- PARAMETERS FOR INTERACTION ---------------'
if (DM_B_ON == 1):


    mx = float(input('Enter the mass of the dark matter particle (in GeV): '))         # Mass of dark matter particle in GeV
    mb = 0.938                         # Mass of baryonic particle in GeV
    sigma_45 = float(input('Enter the interaction cross-section(in units of 10^-45 m^2): '))   # Interaction cross-section in unit of 10^-45 m^2

else:

    mx = 0.4   #doesn't matter
    mb = 0.938                         # Mass of baryonic particle in GeV
    sigma_45 = 0.0
    
print('--------------------------------------------------------------------')    
'-----------------------------------------------------------'



'''
'------------- PRIMORDIAL MAGNETIC FIELD -------------------'

PMF_ON = float(input('Enter 1 to include PMF, 0 to exlude PMF: '))

if(PMF_ON == 1.0):
    B0 = 3.0e-9                         # magnetic field in gauss
    AD_ON = 1.0                         # write 1 to enable ambipolar diffusion heating
    DT_ON = 1.0  
    
else:
    B0 = 3.0e-9                         # magnetic field in gauss
    AD_ON = 0.0                         # write 1 to enable ambipolar diffusion heating
    DT_ON = 0.0
    
'-----------------------------------------------------------'
'''

'----------------------------'
'''-----HELIUM FRACTION-----'''
'----------------------------'

Y = 0.24
f_He = 0.079







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
'We have assumed a matter dominated universe'

def H(red):
    H_z= H0*(Omega_m*(1.0 + red)**3.0)**0.5
    return H_z







'-------------------------------------------------'
'''FUNCTION TO DETERMINE HUBBLE PARAMETER IN CGS'''
'-------------------------------------------------'
'We have assumed a matter dominated universe'

def H_CGS(red):
    return H0_CGS*(Omega_m*(1.0 + red)**3.0)**0.5
 





   
'---------------------------------------------------------'
'''FUNCTION TO DETERMINE NEUTRAL HYDROGEN NUMBER DENSITY'''
'---------------------------------------------------------'

def nH(red):
    NH = 8.403*Omega_b_h2*(1.0 + red)**3.0              #m^-3
    return NH

def nH_CGS(red):
    return 8.403e-6*Omega_b_h2*(1.0 + red)**3.0          #cm-3







'-----------------------------------------------------------'
'''FUNCTION TO DETERMINE RECOMBINATION COEFFICIENT (alpha)'''
'-----------------------------------------------------------'

def alpha_e(Tg): #SI unit
    a = 4.309
    b = -0.6166
    cp = 0.6703
    d = 0.5300
    F = 1.14             # Fudge factor
    t = Tg/(1.0e4)
    
    alpha = F*1.0e-19*((a*t**b)/(1.0 + cp*t**d))    #m^3 s^-1
    return alpha







'------------------------------------------------------------'
'''FUNCTION TO DETERMINE PHOTOIONIZATION COEFFICIENT (beta)'''  
'------------------------------------------------------------'
 
def beta_e(T_gamma): #SI unit, note that T_gamma has been used to calculate beta as suggested in chluba, 2015, mnras
    beta = alpha_e(T_gamma)*2.4093e21*T_gamma**(1.5)*np.exp(-39420.289/T_gamma)
    return beta







'---------------------------'
'''FUNCTION TO DETERMINE C'''
'---------------------------'

def C1(red,x, Tg): #unit less
    K = 7.1898e-23/(H(red))
    Lambda = 8.22458 
    
    Cr = (1.0 + K*Lambda*(1.0 - x)*nH(red))/(1.0 + K*(Lambda + beta_e(Tg))*(1.0 - x)*nH(red))
    return Cr







'----------------------------------------------------'
'''FUNCTION TO DETERMINE DARK MATTER ENERGY DENSITY'''
'----------------------------------------------------'

def rho_DM(red):  #SI unit
    rho_DM_eng =  (Omega_m-Omega_b)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_DM_eng







'-----------------------------------------------'
'''FUNCTION TO DETERMINE MATTER ENERGY DENSITY'''
'-----------------------------------------------'

def rho_M(red): #SI unit
    rho_M_eng = (Omega_m)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_M_eng    







'-----------------------------------------------'
'''FUNCTION TO DETERMINE BARYON ENERGY DENSITY'''
'-----------------------------------------------'

def rho_B(red): #SI unit
    rho_B_eng = (Omega_b)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_B_eng







'------------------------------'
'''FUNCTION TO DETERMINE u_th'''
'------------------------------'

def u_th(Tg, Tx): #SI unit m/s
    return c*2.936835e-7*(((Tg/mb) + (Tx/mx))**(0.5))







'------------------------------'
'''FUNCTION TO DETERMINE F(r)'''
'------------------------------'
  
def F_r(vel_xb, Tg, Tx): #unit less
    u_therm = u_th(Tg, Tx)
    rv = vel_xb/u_therm
    F = sp1.erf(rv/np.sqrt(2.0)) - np.sqrt(2.0/np.pi)*np.exp((-rv**2.0)/2.0)*rv
    return F







'----------------------------------------'
'''FUNCTION TO DETERMINE F(r)/Vel_xb^2'''
'----------------------------------------'
def Fr_by_velxb2(vel_xb, Tg, Tx): #depends on the unit of Vel_xb^2
    u_therm = u_th(Tg, Tx)
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
def Fr_by_velxb(vel_xb, Tg, Tx): #depends on the unit of Vel_xb^2
    u_therm = u_th(Tg, Tx)
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

def Drag_vxb(vel_xb, red, Tg, Tx): #SI unit
    D_Vxb2 = 1.*2.63e+7*h*(Omega_m**0.5)*((1.0 + red)**0.5)*sigma_45*Fr_by_velxb2(vel_xb, Tg, Tx)/(mb + mx)
    return D_Vxb2








'---------------------------------------'
'''FUNCTION TO DETERMINE Q_b_coupling'''
'---------------------------------------'

def Q_b_coupling(Tx, Tg, red, vel_xb): #SI unit
    u_therm = u_th(Tg, Tx)
    rv = vel_xb/u_therm
    Q_b1 = (2.0/3.0)*2.10e+7*(Omega_m - Omega_b)*(h**2.0)*((1.0 + red)**0.5)*sigma_45*np.exp(-(rv**2)/2.0)*(mb/((mb + mx)**2.0))*(Tx - Tg)/(h*(Omega_m**0.5)*u_therm**3)/1.
    return Q_b1






'-----------------------------'
'''FUNCTION TO DETERMINE Q_b_drag'''
'-----------------------------'

def Q_b_drag(Tx, Tg, red, vel_xb): #SI unit
    Q_b_d = 1.*2.26e+3*(Omega_m - Omega_b)*(h**2.)*((1.0+red)**0.5)*sigma_45*(mb*mx/((mb + mx)**2.0))*Fr_by_velxb(vel_xb, Tg, Tx)/(h*Omega_m**0.5) 
    return Q_b_d







'-----------------------------'
'''FUNCTION TO DETERMINE Q_x'''
'-----------------------------'

def Q_x_coupling(Tx, Tg, red, vel_xb): #SI unit
    u_therm = u_th(Tg, Tx)
    rv = vel_xb/u_therm
    Q_x1 = (2.0/3.0)*2.10e+7*Omega_b*(h**2.0)*((1.0 + red)**0.5)*sigma_45*np.exp(-(rv**2)/2.0)*(mx/((mb + mx)**2.0))*(Tg - Tx)/(h*(Omega_m**0.5)*u_therm**3)/1.
    return Q_x1







'----------------------------------'
'''FUNCTION TO DETERMINE Q_x_drag'''
'----------------------------------'

def Q_x_drag(Tx, Tg, red, vel_xb): #SI unit
    Q_x_d = 1.*2.26e+3*Omega_b*h**2*(1.0+red)**0.5*sigma_45*(mb*mx/((mb + mx)**2.0))*Fr_by_velxb(vel_xb, Tg, Tx)/(h*Omega_m**0.5)  
    return Q_x_d







'''----------------------------------------------------------------------------------
FUNCTION QUANTIFIES HEATING DUE TO X-RAY PHOTONS . A POLYNOMIAL FIT TO DATA PROVIDED BY RAGHU(PAGE: 72-73 NOTEBOOK-2) (REF: EQ. 16, PRITCHARD & FURLANETTO, 2006, ASTRO-PH/0607234V2). HEATING FRACTION FROM SHULL AND VAN STEENBERG (1985).
----------------------------------------------------------------------------------'''

def epsilonX(red, x_e):

    #The following fitting formula for heating fraction is taken from Shull and van Steenberg (1985)
    C_heat = 0.9971
    a_heat = 0.2663
    b_heat = 1.3163

    f_heat = C_heat*(1.0 - (1.0 - x_e**a_heat)**b_heat)



    if red > 40.0:
        return 0.0
        
    else:      

        if(source_model == 1): 

            '------- STARBURST GALAXIES (alpha_s = 1.5) --------'
            
		    # A 10th order polynomial is fitted
            
            a0 =  3.0720147690980284e-31      #unit erg/cm^3
            a1 =  -1.6252019812537561e-31
            a2 =  3.6601840659869137e-32
            a3 =  -4.580372420587383e-33
            a4 =  3.5614293381552163e-34
            a5 =  -1.8169143714616605e-35
            a6 =  6.211273713561144e-37
            a7 =  -1.413859771632769e-38
            a8 =  2.0605110229809655e-40
            a9 =  -1.7422341051922115e-42
            a10 =  6.50775481460825e-45
            


        elif(source_model == 2):

            '------- SUPERNOVAE REMNANTS (alpha_s = 1.0) -------'
		
            a0 =  1.3996490397133588e-31      #unit erg/cm^3
            a1 =  -7.378976949697889e-32
            a2 =  1.6571816532400122e-32
            a3 =  -2.069765961650145e-33
            a4 =  1.6071870781337591e-34
            a5 =  -8.191711476559565e-36
            a6 =  2.7985561840414925e-37
            a7 =  -6.367198899518043e-39
            a8 =  9.275888999043167e-41
            a9 =  -7.840764285978355e-43
            a10 =  2.928030972801895e-45


        else:
	
            '------- MINI-QUASARS (alpha_s = 0.5) -------------'
		
            a0 =  5.1293501333535e-32         #unit erg/cm^3
            a1 =  -2.690301413813602e-32
            a2 =  6.020834587151954e-33
            a3 =  -7.506068438235321e-34
            a4 =  5.824551042466467e-35
            a5 =  -2.9689724366878356e-36
            a6 =  1.0148968244730907e-37
            a7 =  -2.3112252188476175e-39
            a8 =  3.3709898359548395e-41
            a9 =  -2.853234621021188e-43
            a10 =  1.0670323740767927e-45
            

        
        epsilonXX = a0 + a1*red + a2*red**2 + a3*red**3 + a4*red**4 + a5*red**5 + a6*red**6 + a7*red**7 + a8*red**8 + a9*red**9 + a10*red**10
        epsilonXX = -2.0*f_star*f_heat*epsilonXX/(3.0*(1.0 + red)*H_CGS(red)*kB_CGS*nH_CGS(red)) 
                                                                                         
        
        return (epsilonXX)







'''----------------------------------------------------------------------------------
FUNCTION QUANTIFIES ionization rate DUE TO X-RAY PHOTONS . PAGE: 72-73 NOTEBOOK-2) (REF: EQ. 16 and subsequent texts, PRITCHARD & FURLANETTO, 2006, ASTRO-PH/0607234V2). IONIZING ENERGY FRACTION FROM SHULL AND VAN STEENBERG (1985).
----------------------------------------------------------------------------------'''

def ionizationrateX(red, x_e):

    Eth = 13.6*1.602e-12 #ionization potential in erg for Hydrogen

    #The following fitting formula for heating fraction is taken from Shull and van Steenberg (1985)
    C_heat = 0.9971
    a_heat = 0.2663
    b_heat = 1.3163

    f_heat = C_heat*(1.0 - (1.0 - x_e**a_heat)**b_heat)

    #The following fitting formula for Hydrogen ionization fraction is taken from Shull and van Steenberg (1985)
    C_ion = 0.3908
    a_ion = 0.4092
    b_ion = 1.7592

    f_ion = C_ion*((1.0 - x_e**a_ion)**b_ion)


    
    if red > 40.0:
        return 0.0
        
    else:        
    
        if(source_model == 1): 

            '------- STARBURST GALAXIES (alpha_s = 1.5) --------'
		
            # A 10th order polynomial is fitted
            
            a0 =  3.0720147690980284e-31      #unit erg/cm^3
            a1 =  -1.6252019812537561e-31
            a2 =  3.6601840659869137e-32
            a3 =  -4.580372420587383e-33
            a4 =  3.5614293381552163e-34
            a5 =  -1.8169143714616605e-35
            a6 =  6.211273713561144e-37
            a7 =  -1.413859771632769e-38
            a8 =  2.0605110229809655e-40
            a9 =  -1.7422341051922115e-42
            a10 =  6.50775481460825e-45

        elif(source_model == 2):

            '------- SUPERNOVAE REMNANTS (alpha_s = 1.0) -------'
		
            a0 =  1.3996490397133588e-31      #unit erg/cm^3
            a1 =  -7.378976949697889e-32
            a2 =  1.6571816532400122e-32
            a3 =  -2.069765961650145e-33
            a4 =  1.6071870781337591e-34
            a5 =  -8.191711476559565e-36
            a6 =  2.7985561840414925e-37
            a7 =  -6.367198899518043e-39
            a8 =  9.275888999043167e-41
            a9 =  -7.840764285978355e-43
            a10 =  2.928030972801895e-45


        else:
	
            '------- MINI-QUASARS (alpha_s = 0.5) -------------'
		
            a0 =  5.1293501333535e-32         #unit erg/cm^3
            a1 =  -2.690301413813602e-32
            a2 =  6.020834587151954e-33
            a3 =  -7.506068438235321e-34
            a4 =  5.824551042466467e-35
            a5 =  -2.9689724366878356e-36
            a6 =  1.0148968244730907e-37
            a7 =  -2.3112252188476175e-39
            a8 =  3.3709898359548395e-41
            a9 =  -2.853234621021188e-43
            a10 =  1.0670323740767927e-45


        epsilonXX = a0 + a1*red + a2*red**2 + a3*red**3 + a4*red**4 + a5*red**5 + a6*red**6 + a7*red**7 + a8*red**8 + a9*red**9 + a10*red**10
        epsilonXX = f_heat*epsilonXX
        ionizationrate = -1.0*f_ion*f_star*epsilonXX/(f_heat*(1.0 + red)*H_CGS(red)*nH_CGS(red)*Eth)   

        
        return (ionizationrate)

'''----------------------------------------------------------------------------------
FUNCTION QUANTIFIES heating rate DUE TO Ly-Alpha PHOTONS. (REF: EQ. 16 and subsequent texts, PRITCHARD & FURLANETTO, 2006, ASTRO-PH/0607234V2). EXCITATION ENERGY FRACTION FROM SHULL AND VAN STEENBERG (1985).
----------------------------------------------------------------------------------'''

def epsilonX_alpha(red, x_e):


    #The following fitting formula for heating fraction is taken from Shull and van Steenberg (1985)
    C_heat = 0.9971
    a_heat = 0.2663
    b_heat = 1.3163

    f_heat = C_heat*(1.0 - (1.0 - x_e**a_heat)**b_heat)


    #The following fitting formula for Hydrogen excitation fraction is taken from Shull and van Steenberg (1985)
    C_exc = 0.4766
    a_exc = 0.2735
    b_exc = 1.5221

    f_exc = C_exc*((1.0 - x_e**a_exc)**b_exc)
    
    if red > 40.0:
        return 0.0
        
    else:      

        if(source_model == 1): 

            '------- STARBURST GALAXIES (alpha_s = 1.5) --------'
		
            # A 10th order polynomial is fitted
            
            a0 =  3.0720147690980284e-31      #unit erg/cm^3
            a1 =  -1.6252019812537561e-31
            a2 =  3.6601840659869137e-32
            a3 =  -4.580372420587383e-33
            a4 =  3.5614293381552163e-34
            a5 =  -1.8169143714616605e-35
            a6 =  6.211273713561144e-37
            a7 =  -1.413859771632769e-38
            a8 =  2.0605110229809655e-40
            a9 =  -1.7422341051922115e-42
            a10 =  6.50775481460825e-45

        elif(source_model == 2):

            '------- SUPERNOVAE REMNANTS (alpha_s = 1.0) -------'
		
            a0 =  1.3996490397133588e-31      #unit erg/cm^3
            a1 =  -7.378976949697889e-32
            a2 =  1.6571816532400122e-32
            a3 =  -2.069765961650145e-33
            a4 =  1.6071870781337591e-34
            a5 =  -8.191711476559565e-36
            a6 =  2.7985561840414925e-37
            a7 =  -6.367198899518043e-39
            a8 =  9.275888999043167e-41
            a9 =  -7.840764285978355e-43
            a10 =  2.928030972801895e-45


        else:
	
            '------- MINI-QUASARS (alpha_s = 0.5) -------------'
		
            a0 =  5.1293501333535e-32         #unit erg/cm^3
            a1 =  -2.690301413813602e-32
            a2 =  6.020834587151954e-33
            a3 =  -7.506068438235321e-34
            a4 =  5.824551042466467e-35
            a5 =  -2.9689724366878356e-36
            a6 =  1.0148968244730907e-37
            a7 =  -2.3112252188476175e-39
            a8 =  3.3709898359548395e-41
            a9 =  -2.853234621021188e-43
            a10 =  1.0670323740767927e-45
        
        epsilonXX = a0 + a1*red + a2*red**2 + a3*red**3 + a4*red**4 + a5*red**5 + a6*red**6 + a7*red**7 + a8*red**8 + a9*red**9 + a10*red**10
        epsilonXX_heat = f_heat*epsilonXX
        
        p_alpha = 0.79
        epsilonXX_alpha = epsilonXX_heat*(f_exc/f_heat)*p_alpha
        
        epsilonXX_alpha = 2.0*f_star*epsilonXX_alpha/(3.0*(1.0 + red)*H_CGS(red)*kB_CGS*nH_CGS(red)) 

        return (epsilonXX_alpha)
        
        


'----------------------------------------------------------------------------------'
'''FUNCTION TO DETERMINE GAS KINETIC TEMPERATURE (Tg) AND IONIZATION FRACTION (x)'''
'----------------------------------------------------------------------------------'

def func(r,red):
    Tg = r[0]
    x = r[1]
    Tx = r[2]
    V_xb = r[3]
    
    f_Tg = ((2.0*Tg)/(1.0 + red)) - ((2.70877e-20*(T_CMB(red) - Tg)*(1.0 + red)**(1.5)*(x/(1.0 + f_He + x)))/(H0*np.sqrt(Omega_m))) - DM_B_ON*Q_b_coupling(Tx, Tg, red, V_xb) - DM_B_ON*Q_b_drag(Tx, Tg, red, V_xb)/1.0 + Xray_ON*epsilonX(red, x) #- Ly_alpha_ON*epsilonX_alpha(red,x)

    f_x =  (C1(red,x,Tg)*(alpha_e(Tg)*x**2*nH(red) - beta_e(T_CMB(red))*(1.0 - x)*np.exp(-118260.87/T_CMB(red))))/(H(red)*(1.0 + red)) + Xray_ON*ionizationrateX(red, x)

    #f_x =  (C1(red,x,Tg)*(alpha_e(Tg)*x**2*nH(red) - beta_e(T_CMB(red))*(1.0 - x)*np.exp(-118260.87/Tg)))/(H(red)*(1.0 + red)) 

    f_Tx = ((2.0*Tx)/(1.0 + red)) - Q_x_coupling(Tx, Tg, red, V_xb) - Q_x_drag(Tx, Tg, red, V_xb)/1.0

    f_Vxb = (V_xb/(1.0 + red)) + Drag_vxb(V_xb, red, Tg, Tx)
    
    return np.array([f_Tg, f_x, f_Tx, f_Vxb], float)






'-------------------------------------------------------------------------------------------------'
'''FUNCTION TO DETERMINE THE VALUE OF PROBABILITY DISTRIBUTION FOR EACH INTIAL RELATIVE VELOCITY'''
'-------------------------------------------------------------------------------------------------'

def prob_func(v_xb):
    return 4*np.pi*v_xb**2*(np.exp((-3*v_xb**2)/(2*Vrms**2)))*(1.0/((2.0/3)*np.pi*Vrms**2)**(1.5))
 







'''----------------------------------------------------------------------------'''
'''--------------------------------EVOLUTION STARTS----------------------------'''
'''----------------------------------------------------------------------------'''


'''---------------REDSHIFT RANGE---------------------------'''
zf = 10.0                       # Final redshift
zi = 1010.0                     # Initial redshift
del_z = -0.01                   # Step-size
z = np.arange(zi, zf, del_z)
'''--------------------------------------------------------'''




Vxb = np.arange(2.0e-6*c, 5.0e-4*c, 2.0e-5*c)       # Range of initial relative velocties
Vrms = 1.0e-4*c                                     # RMS velocity
T_21_vxb = []                                        # Array to store T_21 data for each initial relative velocities (array within and array)
P_Vxb = []                                          # the value of the probabilty distribution for each initial velocity
x_frac_vxb = []                                     # Ionization fraction for each initial relative velocities
T_gas_vxb = []                                      # Gas temperature for each initial relative velocities
T_dark_vxb = []                                     # Dark matter temperature for each initial relative velocites
T_spin_vxb = []
x_a_vxb = []

T_gamma = np.zeros(len(z), float)
K_HH = np.zeros(len(z), float)
nHI = np.zeros(len(z), float)
C_10 = np.zeros(len(z), float)
x_c = np.zeros(len(z), float)
J_a_crit = np.zeros(len(z), float)
J_a_X = np.zeros(len(z), float)
x_a = np.zeros(len(z), float)
T_spin = np.zeros(len(z), float)


#This loop evolves for each intitial relative velocity
for vxb in Vxb:
        

    '----------------------------------------------------'
    '''---------------INITIAL CONDITIONS---------------'''
    '---- (Tg = T_CMB at z = 1010 and x = 0.05497) ----'''
    '-----(Tx = 0.0 at z = 1010 and V_xb_0 = 1e-4*c)---'''
    '----------------------------------------------------'
    
    r0 = np.array([T_CMB(zi), 0.05497, 0.0, vxb], float)    # Initial conditions for Tg, x, Tx and V_xb # 0.95919324





    '''------SOLVING THE EQUATIONS--------'''
    
    r = sp.odeint(func, r0, z)                  # Solving the coupled differential equation
    T_gas = r[:,0]                              # Stores Tg as an array in K
    x_points = r[:,1]                           # Stores x as an array 
    T_dark = r[:,2]                             # Stores Tx as an array in K
    Vel_xb = r[:,3]                             # Stores V_xb as an array in m/s
 

    #T_gamma = T_CMB(z)
    
    for i in range(len(z)):
        
    
    
        '----------CMB TEMPERATURE---------'

        T_gamma[i] = T_CMB(z[i])




        
        '-------CALCULATION OF COLLISIONAL COUPLING COEFFICIENT------'
        
        K_HH[i] = 3.1e-17*T_gas[i]**(0.357)*np.exp(-32.0/T_gas[i])     # Using the fitting formula
        nHI[i] = 8.403*Omega_b_h2*(1.0 + z[i])**3.0
        
        C_10[i] = K_HH[i]*nHI[i]
        
        x_c[i] = (T_star*C_10[i])/(A_10*T_CMB(z[i]))     # Collisional coupling coefficient 



        '-------CALCULATION OF LYMAN-ALPHA COUPLING COEFFICIENT------'


        J_a_crit[i] = 5.5678848e-12*(1.0 + z[i])
        J_a_X[i] = (c_CGS/(4*np.pi))*(Ly_alpha_ON*epsilonX_alpha(z[i], x_points[i]))/(hp_CGS*nu_alpha*(1.0/H_CGS(z[i])*nu_alpha))

        S_a = 1.0

        x_a[i] = 10*(S_a*J_a_X[i])/J_a_crit[i]




        '------CALCULATION OF SPIN TEMPERATURE------'

        #T_spin[i] = ((1.0 +  x_c[i])*T_gas[i]*T_CMB(z[i]))/(x_c[i]*T_gamma[i] + T_gas[i])
        
        T_spin[i] = ((1.0 +  x_c[i] + x_a[i])*T_gas[i]*T_CMB(z[i]))/((x_c[i] + x_a[i])*T_CMB(z[i]) + T_gas[i])
    


    
    T_21 = 0.023*((0.15/Omega_m)*((1.0 + z)/10))**(0.5)*(Omega_b*h/0.02)*(1.0 - (T_gamma/T_spin))  #T_21 brightness temperature in K
    P_Vxb.append(prob_func(vxb))
    T_gas_vxb.append(T_gas)
    T_dark_vxb.append(T_dark)
    T_21_vxb.append(T_21)
    x_a_vxb.append(x_a)
    
    x_frac_vxb.append(x_points)
    T_spin_vxb.append(T_spin)



    #Vel_xb1 = Vel_xb/(1.0 + z)
    #UU_TH = u_th(T_gas, T_dark)
    #rr = Vel_xb/UU_TH
    #Fr = F_r(Vel_xb, T_gas, T_dark)



    
'''------The following steps does a statistical average of the T_21 signal over the velocity distribution---'''
    
T_gas_avg = []
T_dark_avg = []
T_b_avg = []
x_frac_avg = []
T_spin_avg = []
x_a_avg = []


for i in range(len(P_Vxb)):

    T_gas_avg.append(T_gas_vxb[i]*P_Vxb[i])
    T_dark_avg.append(T_dark_vxb[i]*P_Vxb[i])
    T_b_avg.append(T_21_vxb[i]*P_Vxb[i])
    x_frac_avg.append(x_frac_vxb[i]*P_Vxb[i])
    T_spin_avg.append(T_spin_vxb[i]*P_Vxb[i])
    x_a_avg.append(x_a_vxb[i]*P_Vxb[i])

T_gas_avg = sum(T_gas_avg)/(1.0*sum(P_Vxb))
T_dark_avg = sum(T_dark_avg)/(1.0*sum(P_Vxb))
T_21_avg = sum(T_b_avg)/(1.0*sum(P_Vxb))
x_frac_avg = sum(x_frac_avg)/(sum(P_Vxb))
T_spin_avg = sum(T_spin_avg)/(1.0*sum(P_Vxb))
x_a_avg = sum(x_a_avg)/(sum(P_Vxb))



Temp_data = np.array([z, T_gamma, T_gas_avg, T_dark_avg, T_21_avg, x_frac_avg])
Temp_data = Temp_data.T


#Percentage_diff_xfrac=(xfrac_standard_17-x_frac_avg[99280])*100.0/xfrac_standard_17


#print (z[99280], T_gas_avg[99280], T_21_avg[99280], Q_b_coupling(T_dark_avg[99280], T_gas_avg[99280], 17.2, 30.), Q_b_drag(T_dark_avg[99280], T_gas_avg[99280], 17.2, 30.),epsilonX(25.))#, Percentage_diff_xfrac)

'''
if(Xray_ON == 0 and Ly_alpha_ON == 0 and DM_B_ON == 0): 
    np.savetxt('Temp_xfracs_all_zero_standard.txt', Temp_data)

elif(Xray_ON == 0 and Ly_alpha_ON == 0 and DM_B_ON == 1):
    np.savetxt('Temp_xfracs_only_DMb_mx_{0}_sigma_{1}.txt'.format(mx, sigma_45), Temp_data)

elif(Xray_ON == 0 and Ly_alpha_ON == 1 and DM_B_ON == 0):
    np.savetxt('Temp_xfracs_only_Ly_alpha_standard.txt', Temp_data)

elif(Xray_ON == 0 and Ly_alpha_ON == 1 and DM_B_ON == 1):    
    np.savetxt('Temp_xfracs_Ly_alpha_DMb_mx_{0}_sigma_{1}.txt'.format(mx, sigma_45), Temp_data)

elif(Xray_ON == 1 and Ly_alpha_ON == 0 and DM_B_ON == 0):
    np.savetxt('Temp_xfracs_only_Xray_standard_model_{0}.txt'.format(source_model), Temp_data)

elif(Xray_ON == 1 and Ly_alpha_ON == 0 and DM_B_ON == 1):
    np.savetxt('Temp_xfracs_Xray_DMb_mx_{0}_sigma_{1}_model_{2}.txt'.format(mx, sigma_45, source_model), Temp_data)
    
elif(Xray_ON == 1 and Ly_alpha_ON == 1 and DM_B_ON == 0):   
    np.savetxt('Temp_xfracs_Xray_Ly_alpha_standard_model_{0}.txt'.format(source_model), Temp_data)
    
else: 
    np.savetxt('Temp_xfracs_Xray_Ly_alpha_DMb_mx_{0}_sigma_{1}_model_{2}.txt'.format(mx, sigma_45, source_model), Temp_data)
'''    
    
    


print('Redshift = ', z[99280])  # z = 17.2
print('T_21 (K) = ', T_21_avg[99280])
print('Ionization fraction = ', x_frac_avg[99280])
print('Gas Temp (K) = ', T_gas_avg[99280])



'''--------------------  EDGES DATA -------------------------------------'''

EDGES_REDSHIFT, EDGES_T21 = np.loadtxt('EDGES_DATA_T21.txt', unpack = True)

EDGES_REDSHIFT = EDGES_REDSHIFT[12:128]
EDGES_T21 = EDGES_T21[12:128]



'''
plt.figure(1)
plt.loglog(z, T_gamma, color = 'k', linestyle = '--', label = 'CMB')
plt.loglog(z, T_gas_avg, color = 'b', label = 'IGM')
#plt.loglog(z, T_dark_avg, color = 'r', label = 'DM')
plt.loglog(z, T_spin_avg, color = 'r', label = 'Spin')
plt.xlabel('z')
plt.ylabel('Temperature (in K)')
plt.title('Temperature evolution')
plt.legend() 
plt.savefig('Temp_evolve.png', dpi = 300)'''

'''
plt.figure(2)
plt.loglog(z, x_frac_avg)
plt.xlabel('z')
plt.ylabel(r'Ionization fraction $x$')
plt.title('Ionization fraction evolution')
'''


plt.figure(3)
plt.plot(z, T_21_avg, label = r'Our model $f_*$ = {0}'.format(f_star*0.01))
#plt.plot(EDGES_REDSHIFT, EDGES_T21, label = 'EDGES')
plt.xlim(10.0, 300.0)
#plt.axvline(16.8, color = 'k', linestyle = '--')
#plt.axvline(16.2, color = 'g', linestyle = '--')
#plt.axvline(17.2, color = 'r', linestyle = '--')
plt.xlabel('z')
plt.ylabel(r'$T_{21}$ (in K)')
plt.legend()
plt.title(r'$f_*$ = {0} , Green = 16.2 & Black = 16.8 & Red = 17.2'.format(f_star*0.01))

#plt.savefig('T_21_Rise_f_star0p01_sigma100', dpi = 150)
'''
plt.figure(4)
plt.loglog(z, x_a_avg)
plt.xlabel('z')
'''
plt.show()
