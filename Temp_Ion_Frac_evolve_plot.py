
import matplotlib.pylab as plt
import matplotlib.pyplot as pyplt
import numpy as np
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mticker


source_model = 1

if (source_model == 1):


    data_100 = np.loadtxt('Temp_xfracs_mx0.4GeV_sigma100_STARBURST.txt', float)
    data_0p2 = np.loadtxt('Temp_xfracs_mx0.4GeV_sigma0.2_STARBURST.txt', float)
    data_standard = np.loadtxt('Temp_xfracs_standard_STARBURST.txt', float)


elif(source_model == 2):

    data_100 = np.loadtxt('Temp_xfracs_mx0.4GeV_sigma100_SNR.txt', float)
    data_0p2 = np.loadtxt('Temp_xfracs_mx0.4GeV_sigma0.2_SNR.txt', float)
    data_standard = np.loadtxt('Temp_xfracs_standard_SNR.txt', float)



else:

    data_100 = np.loadtxt('Temp_xfracs_mx0.4GeV_sigma100_MINIQUASARS.txt', float)
    data_0p2 = np.loadtxt('Temp_xfracs_mx0.4GeV_sigma0.2_MINIQUASARS.txt', float)
    data_standard = np.loadtxt('Temp_xfracs_standard_MINIQUASARS.txt', float)





z = data_100[:,0]
T_gamma = data_100[:,1]
T_gas_avg_100 = data_100[:,2]
T_dark_avg_100 = data_100[:,3]
T_21_avg_100 = data_100[:,4]
x_points_avg_100 = data_100[:,5]


T_gas_avg_0p2 = data_0p2[:,2]
T_dark_avg_0p2 = data_0p2[:,3]
T_21_avg_0p2 = data_0p2[:,4]
x_points_avg_0p2 = data_0p2[:,5]

T_gas_standard = data_standard[:,2]
x_points_standard = data_standard[:,5]




'''----------------------------------------------------------------'''
'''---------------------PLOT TEMPERATURE---------------------------'''
'''----------------------------------------------------------------'''


plt.loglog(z, T_gamma, 'k:')#, label = r'$T_{CMB}$')
plt.loglog(z, T_gas_standard, 'k', label = r'$T_g$ (standard)')
plt.loglog(z, T_gas_avg_100, 'r--', label = r'$m_{\chi}$ = 0.4 GeV, $\sigma_{45}$ = 100 ')
plt.loglog(z, T_dark_avg_100, 'r--')
plt.loglog(z, T_gas_avg_0p2, 'b-.', label = r'$m_{\chi}$ = 0.4 GeV, $\sigma_{45}$ = 0.2 ')
plt.loglog(z, T_dark_avg_0p2, 'b-.')



plt.xlim(1e1, 1e3)
plt.ylim(1e0, 5e3)

plt.axvspan(15, 22, linewidth = 0.5, alpha = 0.7, color = 'salmon')

plt.xlabel('z')
plt.ylabel('Temperature (K)')
plt.legend(loc='best')

if(source_model == 1):

    #plt.title('Starburst Galaxies')     # Use this when plotting for STARBURST GALAXIES
    #plt.savefig('Temperature_evolution_STARBURST.png', dpi = 400)

elif(source_model == 2):

    plt.title('Supernovae Remnants')    # Use this when plotting for STARBURST GALAXIES
    #plt.savefig('Temperature_evolution_SNR.png', dpi = 400)

else:

    plt.title('Miniquasars')            # Use this when plotting for MINIQUASARS
    #plt.savefig('Temperature_evolution_MINIQUASARS.png', dpi = 400)






'''------------------------------------------------------------------------'''
'''---------------------PLOT IONIZATION FRACTION---------------------------'''
'''------------------------------------------------------------------------'''


percentage_diff_x_mx_sigma_100 = abs(x_points_avg_100 - x_points_standard)*100.0/x_points_standard
percentage_diff_x_mx_sigma_0p2 = abs(x_points_avg_0p2 - x_points_standard)*100.0/x_points_standard

grid = plt.GridSpec(2, 1, height_ratios=[3, 2])
 
#pyplt.tight_layout()
fig = pyplt.figure(figsize=(5,5))
 
ax0 = plt.subplot(grid[0])
line0, = plt.loglog(z, x_points_standard/1.0e-3, 'k-')
line1, = plt.loglog(z, x_points_avg_100/1.0e-3, 'r--')
line2, = plt.loglog(z, x_points_avg_0p2/1.0e-3, 'b-.')
plt.yscale('linear')
plt.xlim(10.0, 200.0)
plt.setp(ax0.get_xticklabels(), visible=False)
 
plt.ylim(0.05e0, 1.0e0)
plt.ylabel(r'$x$', size=12)
plt.title(r'x $10^{-3}$', loc = 'left', size = 10)  



 
ax1 = plt.subplot(grid[1], sharex = ax0)
plt.loglog(z, percentage_diff_x_mx_sigma_100, 'r--')
plt.loglog(z, percentage_diff_x_mx_sigma_0p2, 'b-.')
plt.yscale('linear')
#plt.setp(ax0.get_xticklabels(), visible=False)
 
plt.xlim(1e1, 200.0)
#ax1.set_xticks([10,20,40,70,100,300])
 
#ax1.set_yticks([2,4,6,8,10,12,14,16,18,20])
 
 
plt.ylim(-2.0, 30.0)
 
 
plt.xlabel('z', size=12)
plt.ylabel(r'% change in $x$', size=12)
 
# put lened on first subplot
ax0.legend((line0, line1, line2), ('Standard',r"$m_{\chi}=0.4 \, {\rm GeV} \, , \sigma_{45} = 100$", r"$m_{\chi}=0.4 \, {\rm GeV} \, , \sigma_{45} = 0.2$"), loc='upper right')
 
plt.subplots_adjust(hspace=.0)
#ticks = [0, 10, 20,30,40] 
#plt.set_yticks(ticks)
plt.text(172.0, -7.0, '200')



if(source_model == 1):

    #plt.suptitle('Starburst Galaxies')     # Use this when plotting for STARBURST GALAXIES
    #plt.savefig('Ionization_fraction_evolution_STARBURST.png', dpi = 400)

elif(source_model == 2):

    plt.suptitle('Supernovae Remnants')    # Use this when plotting for STARBURST GALAXIES
    #plt.savefig('Ionization_fraction_evolution_SNR.png', dpi = 400)

else:

    plt.suptitle('Miniquasars')            # Use this when plotting for MINIQUASARS
    #plt.savefig('Ionization_fraction_evolution_MINIQUASARS.png', dpi = 400)


plt.show()




