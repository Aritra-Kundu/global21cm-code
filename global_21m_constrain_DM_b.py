
#This code just solves the coupled equations for a range of mass and cross section. 
#It generates the plots of del_x vs Tb and the constraint on mass and cross-section
#(with averaging over velocity) as given in Munoz's (2015) paper. 

from matplotlib import ticker, cm
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







'----------------------------------------------------------'

E1s = 2.176e-11                     # in erg
Xray_ON = 1.0                       # write '1' to enable decaying turbulence heating, '0' to disable it

source_model = 1                    # 1 = STARBURST ! 2 = SNR ! 3 = MINIQUASARS !
f_star = 1.0         #in units of 0.01      # The data of Xray heating by Raghu is for f_star = 0.01. Change the value accordingly







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
    NH = 8.403*Omega_b_h2*(1.0 + red)**3.0               # m^-3
    return NH

def nH_CGS(red):
    return 8.403e-6*Omega_b_h2*(1.0 + red)**3.0          # cm^-3








'-----------------------------------------------------------'
'''FUNCTION TO DETERMINE RECOMBINATION COEFFICIENT (alpha)'''
'-----------------------------------------------------------'

def alpha_e(Tg):     # SI unit
    a = 4.309
    b = -0.6166
    cp = 0.6703
    d = 0.5300
    F = 1.14             # Fudge factor
    t = Tg/(1.0e4)
    
    alpha = F*1.0e-19*((a*t**b)/(1.0 + cp*t**d))   # m^3 s^-1
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

def C1(red,x, Tg):   # unit less
    K = 7.1898e-23/(H(red))
    Lambda = 8.22458 
    
    Cr = (1.0 + K*Lambda*(1.0 - x)*nH(red))/(1.0 + K*(Lambda + beta_e(Tg))*(1.0 - x)*nH(red))
    return Cr







'---------------------------------------------'
'''FUNCTION TO DETERMINE DARK MATTER DENSITY'''
'---------------------------------------------'

def rho_DM(red):   # SI unit
    rho_DM_eng =  (Omega_m-Omega_b)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_DM_eng







'----------------------------------------'
'''FUNCTION TO DETERMINE MATTER DENSITY'''
'----------------------------------------'

def rho_M(red):   # SI unit
    rho_M_eng = (Omega_m)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_M_eng    







'----------------------------------------'
'''FUNCTION TO DETERMINE BARYON DENSITY'''
'----------------------------------------'

def rho_B(red):   # SI unit
    rho_B_eng = (Omega_b)*((3*H0**2*c**2)/(8*np.pi*G))*(1.0 + red)**3
    return rho_B_eng







'------------------------------'
'''FUNCTION TO DETERMINE u_th'''
'------------------------------'

def u_th(Tg, Tx, dm):   #SI unit m/s
    mx = dm
    return c*2.936835e-7*(((Tg/mb) + (Tx/mx))**(0.5))







'------------------------------'
'''FUNCTION TO DETERMINE F(r)'''
'------------------------------'
  
def F_r(vel_xb, Tg, Tx, dm):   #unit less
    u_therm = u_th(Tg, Tx, dm)
    rv = vel_xb/u_therm
    F = sp1.erf(rv/np.sqrt(2.0)) - np.sqrt(2.0/np.pi)*np.exp((-rv**2.0)/2.0)*rv
    return F







'----------------------------------------'
'''FUNCTION TO DETERMINE F(r)/Vel_xb^2'''
'----------------------------------------'

def Fr_by_velxb2(vel_xb, Tg, Tx, dm):   #depends on the unit of Vel_xb^2
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

def Fr_by_velxb(vel_xb, Tg, Tx, dm):    #depends on the unit of Vel_xb^2
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

def Drag_vxb(vel_xb, red, Tg, Tx, sigma_45, dm):  #SI unit
    mx = dm
    D_Vxb2 = 1.*2.63e+7*h*(Omega_m**0.5)*((1.0 + red)**0.5)*sigma_45*Fr_by_velxb2(vel_xb, Tg, Tx, dm)/(mb + mx)
    return D_Vxb2







'---------------------------------------'
'''FUNCTION TO DETERMINE Q_b_coupling'''
'---------------------------------------'

def Q_b_coupling(Tx, Tg, red, vel_xb, sigma_45, dm):   #SI unit
    mx = dm
    u_therm = u_th(Tg, Tx, dm)
    rv = vel_xb/u_therm
    Q_b1 = (2.0/3.0)*2.10e+7*(Omega_m - Omega_b)*(h**2.0)*((1.0 + red)**0.5)*sigma_45*np.exp(-(rv**2)/2.0)*(mb/((mb + mx)**2.0))*(Tx - Tg)/(h*(Omega_m**0.5)*u_therm**3)/1.0
    return Q_b1








'-----------------------------'
'''FUNCTION TO DETERMINE Q_b_drag'''
'-----------------------------'

def Q_b_drag(Tx, Tg, red, vel_xb, sigma_45, dm):    #SI unit
    mx = dm
    Q_b_d = 1.*2.26e+3*(Omega_m - Omega_b)*h**2*(1.0+red)**0.5*sigma_45*(mb*mx/((mb + mx)**2.0))*Fr_by_velxb(vel_xb, Tg, Tx, dm)/(h*Omega_m**0.5) 
    return Q_b_d







'-----------------------------'
'''FUNCTION TO DETERMINE Q_x_coupling'''
'-----------------------------'

def Q_x_coupling(Tx, Tg, red, vel_xb, sigma_45, dm):   #SI unit
    mx = dm
    u_therm = u_th(Tg, Tx, dm)
    rv = vel_xb/u_therm
    Q_x1 = (2.0/3.0)*2.10e+7*Omega_b*(h**2.0)*((1.0 + red)**0.5)*sigma_45*np.exp(-(rv**2)/2.0)*(mx/((mb + mx)**2.0))*(Tg - Tx)/(h*(Omega_m**0.5)*u_therm**3)/1.0
    return Q_x1








'-----------------------------'
'''FUNCTION TO DETERMINE Q_x_drag'''
'-----------------------------'

def Q_x_drag(Tx, Tg, red, vel_xb, sigma_45, dm):    #SI unit
    mx = dm
    Q_x_d = 1.*2.26e+3*Omega_b*h**2*(1.0+red)**0.5*sigma_45*(mb*mx/((mb + mx)**2.0))*Fr_by_velxb(vel_xb, Tg, Tx, dm)/(h*Omega_m**0.5)  
    return Q_x_d







'''----------------------------------------------------------------------------------
FUNCTION QUNTIFIES HEATING DUE TO X-RAY PHOTONS . A POLYNOMIAL FIT TO DATA PROVIDED BY RAGHU(PAGE: 72-73 NITEBOOK-2) (REF: EQ. 16, PRITCHARD & FURLANETTO, 2006, ASTRO-PH/0607234V2). HEATING FRACTION FROM SHULL AND VAN STEENBERG (1985).
----------------------------------------------------------------------------------'''

def epsilonX(red, x_e):


    #The following fitting formula for heating fraction is taken from Shull and van Steenberg (1985)
    C_heat = 0.9971
    a_heat = 0.2663
    b_heat = 1.3163

    f_heat = C_heat*(1.0 - (1.0 - x_e**a_heat)**b_heat)

    if red > 25.0:
        return 0.0
        
    else:      

        if(source_model == 1): 

            '------- STARBURST GALAXIES (alpha_s = 1.5) --------'
		
            a0 = -9.09011231e-32       #(5 order)  #unit erg/cm^3
            a1 = 2.73837019e-32        #(5 order)
            a2 = -2.92872944e-33       #(5 order)
            a3 = 1.47016982e-34        #(upto a5)
            a4 = -3.54781625e-36       #(upto a5)
            a5 = 3.33317827e-38        #(upto a5)

        elif(source_model == 2):

            '------- SUPERNOVAE REMNANTS (alpha_s = 1.0) -------'
		
            a0 = -4.02235280e-32        #(5 order)  #unit erg/cm^3
            a1 = 1.20940822e-32         #(5 order)
            a2 = -1.29105369e-33        #(5 order)
            a3 = 6.46771200e-35         #(upto a5)
            a4 = -1.55724782e-36       #(upto a5)
            a5 = 1.45926920e-38         #(upto a5)

        else:
	
            '------- MINI-QUASARS (alpha_s = 0.5) -------------'
		
            a0 = -1.39465265e-32       #(5 order)  #unit erg/cm^3
            a1 = 4.19350302e-33        #(5 order)
            a2 = -4.47281812e-34       #(5 order)
            a3 = 2.23764323e-35        #(upto a5)
            a4 = -5.37784562e-37       #(upto a5)
            a5 = 5.02811663e-39        #(upto a5)

        
        epsilonXX = a0 + a1*red + a2*red**2 + a3*red**3 + a4*red**4 + a5*red**5
        epsilonXX = -2.0*f_star*f_heat*epsilonXX/(3.0*(1.0 + red)*H_CGS(red)*kB_CGS*nH_CGS(red)) 
                                                                                         
        
        return (epsilonXX)







'''----------------------------------------------------------------------------------
FUNCTION QUNTIFIES ionization rate DUE TO X-RAY PHOTONS . PAGE: 72-73 NITEBOOK-2) (REF: EQ. 16 and subsequent texts, PRITCHARD & FURLANETTO, 2006, ASTRO-PH/0607234V2). IONIZING ENERGY FRACTION FROM SHULL AND VAN STEENBERG (1985).
----------------------------------------------------------------------------------'''

def ionizationrateX(red, x_e):

    Eth = 13.6*1.602e-12 #ionization potential in erg for Hydrogen
    

    #The following fitting formula for heating fraction is taken from Shull and van Steenberg (1985)
    C_heat = 0.9971
    a_heat = 0.2663
    b_heat = 1.3163

    f_heat = C_heat*(1.0 - (1.0 - x_e**a_heat)**b_heat)

    #The following fitting formula for heating fraction is taken from Shull and van Steenberg (1985)
    C_ion = 0.3908
    a_ion = 0.4092
    b_ion = 1.7592

    f_ion = C_ion*((1.0 - x_e**a_heat)**b_heat)


    
    if red > 25.:
        return 0.0
        
    else:        
    
        if(source_model == 1): 

            '------- STARBURST GALAXIES (alpha_s = 1.5) --------'
		
            a0 = -9.09011231e-32       #(5 order)  #unit erg/cm^3
            a1 = 2.73837019e-32        #(5 order)
            a2 = -2.92872944e-33       #(5 order)
            a3 = 1.47016982e-34        #(upto a5)
            a4 = -3.54781625e-36       #(upto a5)
            a5 = 3.33317827e-38        #(upto a5)

        elif(source_model == 2):

            '------- SUPERNOVAE REMNANTS (alpha_s = 1.0) -------'
		
            a0 = -4.02235280e-32        #(5 order)  #unit erg/cm^3
            a1 = 1.20940822e-32         #(5 order)
            a2 = -1.29105369e-33        #(5 order)
            a3 = 6.46771200e-35         #(upto a5)
            a4 = -1.55724782e-36       #(upto a5)
            a5 = 1.45926920e-38         #(upto a5)

        else:
	
            '------- MINI-QUASARS (alpha_s = 0.5) -------------'
		
            a0 = -1.39465265e-32       #(5 order)  #unit erg/cm^3
            a1 = 4.19350302e-33        #(5 order)
            a2 = -4.47281812e-34       #(5 order)
            a3 = 2.23764323e-35        #(upto a5)
            a4 = -5.37784562e-37       #(upto a5)
            a5 = 5.02811663e-39        #(upto a5)



        epsilonXX = a0 + a1*red + a2*red**2 + a3*red**3 + a4*red**4 + a5*red**5
        epsilonXX = f_heat*epsilonXX
        ionizationrate = -1.0*f_ion*f_star*epsilonXX/(f_heat*(1.0 + red)*H_CGS(red)*nH_CGS(red)*Eth)   

        
        return (ionizationrate)







'----------------------------------------------------------------------------------'
'''FUNCTION TO DETERMINE GAS KINETIC TEMPERATURE (Tg) AND IONIZATION FRACTION (x)'''
'----------------------------------------------------------------------------------'

def func(r,red, cross, dm):
    Tg = r[0]
    x = r[1]
    Tx = r[2]
    V_xb = r[3]
    
    
    f_Tg = ((2.0*Tg)/(1.0 + red)) - ((2.70877e-20*(T_CMB(red) - Tg)*(1.0 + red)**(1.5)*(x/(1.0 + x)))/(H0*np.sqrt(Omega_m))) - Q_b_coupling(Tx, Tg, red, V_xb, cross, dm)-Q_b_drag(Tx, Tg, red, V_xb, cross, dm)/1.0 + Xray_ON*epsilonX(red, x)

    f_x =  (C1(red,x,Tg)*(alpha_e(Tg)*x**2*nH(red) - beta_e(T_CMB(red))*(1.0 - x)*np.exp(-118260.87/Tg)))/(H(red)*(1.0 + red)) + ionizationrateX(red, x)

    #f_x =  (C1(red,x,Tg)*(alpha_e(Tg)*x**2*nH(red) - beta_e(T_CMB(red))*(1.0 - x)*np.exp(-118260.87/Tg)))/(H(red)*(1.0 + red)) 

    f_Tx = ((2.0*Tx)/(1.0 + red)) - Q_x_coupling(Tx, Tg, red, V_xb, cross, dm) -Q_x_drag(Tx, Tg, red, V_xb, cross, dm)/1.0

    f_Vxb = (V_xb/(1.0 + red)) + Drag_vxb(V_xb, red, Tg, Tx, cross, dm)
    
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






'----------CMB TEMPERATURE---------'

T_gamma = T0*(1.0 + z)







Vxb = np.arange(2.0e-6*c, 5.0e-4*c, 2.0e-5*c)            # Range of initial relative velocties
Vrms = 1e-4*c                                            # RMS velocity





N_grid = 100             # Number of grid points

cross_sec_45 = 10**(np.linspace(-2, 3, N_grid))
dark_mass = 10**(np.linspace(-4, 2, N_grid))





sigma_300 = []                # Stores the values of sigma_45 for which T_21 lies within the EDGES bounds
mass_dark_300 = []            # Stores the values of dark matter mass for which T_21 lies within the EDGES bounds






mb = 0.938            # Mass of baryon


'''------------------------------------------------------'''
'''----------------- STANDARD x_frac at 17.2 ------------'''
'''------------------------------------------------------'''

if (source_model == 1):     # STARBURST

    xfrac_standard_17 = 2.108985732768907471e-04

elif (source_model == 2):  # SNR

    xfrac_standard_17 = 2.019655197149658552e-04

else:   # MINIQUASARS

    xfrac_standard_17 = 1.973069325846223577e-04







T_21_density = np.empty([N_grid, N_grid], float)         # Will be used to create density/contour plots of T_21
x_density = np.empty([N_grid, N_grid], float)            # Will be used to create density/contour plots of x



T_bright_avg = []      # Will be used to create the scatter plot ! Stores the T_b for different mass and cross-section
Ion_avg = []           # Will be used to create the scatter plot ! Stores the del_x % for different mass and cross-section







for i in range(0, N_grid, 1):      #Loop for cross_section
    for j in range(0, N_grid, 1):  #Loop for dark matter mass


        
        T_21_xb = []           # Stores T_21 for each parameter set
        P_Vxb = []             
        x_points_xb = []       # Stores x for each parameter set


        
        for vxb in Vxb:



        
            '----------------------------------------------------'
            '''---------------INITIAL CONDITIONS---------------'''
            '---- (Tg = T_CMB at z = 1010 and x = 0.05497)-------'
            '----------------------------------------------------'



            r0 = np.array([T_CMB(zi), 0.05497, 0.0, vxb], float)    # Initial conditions for Tg and x





            '''------SOLVING THE EQUATIONS--------'''

                         
            r = sp.odeint(func, r0, z, args = (cross_sec_45[i], dark_mass[j]))
            T_gas = r[:,0]                        # Stores Tg as an array
            x_points = r[:,1]                     # Stores x as an array
            T_dark = r[:,2]                       # Stores Tx as an array
            Vel_xb = r[:,3]                       # Stores V_xb as an array




            T_21 = 23*((0.15/Omega_m)*((1.0 + z)/10))**(0.5)*(Omega_b*h/0.02)*(1.0 - (T_gamma/T_gas))     # in mK
            P_Vxb.append(prob_func(vxb))
            T_21_xb.append(T_21)
            x_points_xb.append(x_points)



        T_b_avg = []
        x_avg = []



        for o in range(len(P_Vxb)):
            T_b_avg.append(T_21_xb[o]*P_Vxb[o])
            x_avg.append(x_points_xb[o]*P_Vxb[o])




        T_21_avg = sum(T_b_avg)/(1.0*sum(P_Vxb))           # The average T_21 signal for each parameter space (The global signal)
        x_points_avg = sum(x_avg)/(1.0*sum(P_Vxb))         # The average ionization_fraction x for each parameter space




        #T_21_density[i,j] = np.log10(-1*T_21_avg[99280])

        if(-1*T_21_avg[99280] >= 300.0 and -1*T_21_avg[99280] <= 1000.0):   # !!! 99280 = index of z = 17.2 !!!

            sigma_300.append(cross_sec_45[i])
            mass_dark_300.append(dark_mass[j])

            Percentage_diff_xfrac=(xfrac_standard_17-x_points_avg[99280])*100.0/xfrac_standard_17

            x_density[i,j] = abs((xfrac_standard_17-x_points_avg[99280])*100.0/xfrac_standard_17)
            T_21_density[i,j] = np.log10(-1*T_21_avg[99280])

            #print(dark_mass[j],cross_sec_45[i], T_gas[99280],T_21_avg[99280], x_points_avg[99280], Percentage_diff_xfrac)

            T_bright_avg.append(-1*T_21_avg[99280])

            Ion_avg.append(Percentage_diff_xfrac)
        
        else:

            x_density[i,j] = 50.0         # If not with EDGES limit then assign some random value
            T_21_density[i,j] = 50.0      # If not with EDGES limit then assign some random value

        print(i, ',', j)


'''value = 50.0
masked_array = np.ma.masked_where(x_density == 50.0, x_density)
cmap = cm.magma
cmap.set_bad(color = 'white')'''




Temp_x_data = np.array([T_bright_avg, Ion_avg])
Temp_x_data = Temp_x_data.T

contour_data = np.array([dark_mass, cross_sec_45])
contour_data = contour_data.T



if(source_model == 1):

    np.savetxt('Temp_del_x_data_STARBURST.txt', Temp_x_data)
    np.savetxt('Mass_sigma_grid100_contour_STARBURST.txt', contour_data)
    np.savetxt('Ionization_frac_contour_STARBURST.txt', x_density)
    np.savetxt('T_21_contour_STARBURST.txt', T_21_density)

elif(source_model == 2):

    np.savetxt('Temp_del_x_data_SNR.txt', Temp_x_data)
    np.savetxt('Mass_sigma_grid100_contour_SNR.txt', contour_data)
    np.savetxt('Ionization_frac_contour_SNR.txt', x_density)
    np.savetxt('T_21_contour_SNR.txt', T_21_density)


else:

    np.savetxt('Temp_del_x_data_MINIQUASARS.txt', Temp_x_data)
    np.savetxt('Mass_sigma_grid100_contour_MINIQUASARS.txt', contour_data)
    np.savetxt('Ionization_frac_contour_MINIQUASARS.txt', x_density)
    np.savetxt('T_21_contour_MINIQUASARS.txt', T_21_density)




'''#plt.imshow(masked_array, origin = 'lower', cmap = cmap, aspect = 6.0/5.0, interpolation = 'bilinear', vmin = None, vmax = None, extent = [-4, 2, -2, 3])
plt.contourf(np.log10(dark_mass), np.log10(cross_sec_45), x_density, levels = [0, 5, 10, 15, 20, 25, 30, 35])
clb = plt.colorbar()
plt.xlabel(r'$log_{10}(m_{\chi}/{\rm Gev})$')
plt.ylabel(r'$log_{10}(\sigma_{45})$')
#plt.title('Bounds on cross-section and dark matter mass')
clb.set_label(r'$\Delta x$ %')
#plt.savefig('Ionize_contour_3.png', dpi = 400)'''






plt.scatter(T_bright_avg, Ion_avg, c = 'b', s=5, alpha = 1, linewidths=1)
plt.xlabel(r'-$T_b$ (mK)')
plt.ylabel(r'$\Delta$ x %')
plt.xlim(200, 1100)
plt.ylim(-5, 37)
#plt.savefig('Ionize_Tb_80.png', dpi = 400)
plt.show()









