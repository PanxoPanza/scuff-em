# -*- coding: utf-8 -*-
"""
This script compares the plots of scuff-avescatter test for prolate and oblate 
spheroids, comparing with the exact solution

Created on Wed Nov 10 09:43:21 2021

@author: Francisco Ramirez
"""
import numpy as np
import matplotlib.pyplot as plt

spheroid_type = 'oblate' # change name to 'prolate' if needed 
namefile_exact = spheroid_type + "_exact.dat"
namefile_sim =  spheroid_type + ".AVSCAT.EMTPFT"
lnstyle = ['solid','dashed']

sol = []
sol.append(np.loadtxt(namefile_exact))
sol.append(np.loadtxt(namefile_sim))

fig, ax1 = plt.subplots()
ax2 = plt.twinx()
ax2.set_ylim(-1,1)
plt.xscale("log")
ax1.set_xlabel("Wavelength ($\mu$m)")
ax1.set_ylabel(r"Cross sections ($\mu$m$^2$), $\langle C_\mathrm{abs}\rangle$ and $\langle C_\mathrm{sca}\rangle$")
ax2.set_ylabel(r"Asymmetry parameter, $\langle \mu_\mathrm{sca}\rangle$")
plt.show()

# plot exact and simulations
for i in range(len(sol)):
    omega = sol[i][:,1]
    Cabs = sol[i][:,3]
    Csca = sol[i][:,4]
    gcos = sol[i][:,5]/sol[i][:,4]
    ax1.plot(2*np.pi/omega,Cabs,linestyle = lnstyle[i], color='b',
             label = r'$\langle C_\mathrm{abs}\rangle$')
    ax1.plot(2*np.pi/omega,Csca,linestyle = lnstyle[i], color='r',
             label = r'$\langle C_\mathrm{sca}\rangle$')
    ax2.plot(2*np.pi/omega,gcos,linestyle = lnstyle[i], color='k',
             label = r'$\langle \mu_\mathrm{sca}\rangle$')

# Add auxiliary legends
p1, = ax1.plot(-1,-1,'-r')
p2, = ax1.plot(-1,-1,'-b')
p3, = ax1.plot(-1,-1,'-k')
l1 = ax1.legend([p1, p2, p3], [r'$\langle C_\mathrm{abs}\rangle$', 
                           r'$\langle C_\mathrm{sca}\rangle$',
                           r'$\langle \mu_\mathrm{sca}\rangle$'],
                frameon=False)

# second legend
p4, = ax1.plot(-1,-1,'-k')
p5, = ax1.plot(-1,-1,'--k')
l2 = ax1.legend([p4, p5], ['Exact', 'Simulation'],
                frameon=False,
                loc='lower left')

ax1.add_artist(l1)
