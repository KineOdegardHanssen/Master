from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

# Numerical results
# Feed in the spin correlation
# Different default options to make it easier
#infile = open("/home/ubu/Master/Spincorrelation_Results/Spincorrelations_from_simulations/0aa0_fcc20x20x20_beta1_Jsm1_eqsteps10000_mcsteps_inbin_100000_no_of_bins100_spincorrelationfunction.txt", "r")
# 14.03.17
# Made on this computer
#infile = open("/home/ubu/Master/Spincorrelation_Results/Results_test_140317/0aa0_fcc6x6x6_beta0p1_Jxy1_Jyz1_sianDx1_sinDy1_eqsteps10000_mcsteps_inbin_1000_no_of_bins100_spincorrelationfunction.txt", "r")
# Made on the stationary one
#infile = open("/home/ubu/Master/Spincorrelation_Results/Spincorrelations_from_simulations/Js1/fcc6x6x6_Jxym1_Jxz0_Jyz1_sianDx1Dy1Dz0_beta1_eq10000_mc10000_bins100_spincorrelationfunction.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Spincorrelations_from_simulations/Js1/fcc6x6x6_Jxy1_Jxz0_Jyz0_beta10_eq10000_mc1000_bins100_II_spincorrelationfunction.txt", "r")
# Made 15.03.17 on my stationary
#infile = open("/home/ubu/Master/Spincorrelation_Results/Results_test_140317/0aa0_fcc6x6x6_Jxy1_Jyz1_sianDx1_sianDy1_beta1p025_eqsteps10000_mcsteps_inbin_1000_no_of_bins100_spincorrelationfunction.txt", "r")

#infile = open("/home/ubu/Master/Spincorrelation_Results/Spincorrelations_from_simulations/Js1/test2_spincorrelationfunction.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Spincorrelations_from_simulations/0000test_spincorrelationfunction.txt", "r")

#infile = open("/home/ubu/Master/Spincorrelation_Results/Spincorrelations_from_simulations/Realistic_parameters/fcc6x6x6_Jxy1_Jxz0_Jyz1_beta10_eq10000_mc1000_bins100_seed59_spincorralldir_spincorrelationfunctiontot.txt", "r")


# Made on 04.06.17. To test Jensen's parameters
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x8yopen_nnJyz1p04_nnnJy0p67_T1p5K_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc6x6x6_nnJyz1p04_nnnJy0p67_beta100_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_spincorrelationfunctiontot.txt", "r")

# Made on 02.08.17 for further testing
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x8yopen_nnJyz1_nnJxy1_T15K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x8yopen_nnJyz1p33_nnnJy0p67_T1p5K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x8yopen_nnJyz1p04_nnJxy0p3_nnnJy0p67_T1p5K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x8yopen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_T100K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc6x6x6yopen_nnJyz1_nnJxy1_T15K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc6x6x6yopen_nnJyz1_nnJxy1_T15K_eq10000_mc10000_bins100_seed79_latticeseed21_spincorrelationfunctiontot.txt", "r")

#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc6x6x6yopen_nnJyz1_nnJxy1_beta10_eq10000_mc10000_bins100_seed79_latticeseed21_spincorrelationfunctiontot.txt", "r")

#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc6x6x6_periodic_nnJyz1_nnJxy1_T15K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")

#Different length of the dimensions
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x4yopen_nnJyz1p04_nnnJy0p67_T1p5K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")

# Made on 03.08.17 for further testing
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x8yopen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p17_Dy0p91_T14K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")

# Made on 04.08.17 for further testing
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc12x12x12yopen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_T21K_eq10000_mc1000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")

#fcc8x8x8yopen_nnJyz1p04_nnnJy0p67_T18p5K_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_spincorrelationfunctiontot
# fcc6x6x6_nnJyz1p04_nnnJy0p67_beta0p1_eq10000_mc10000_bins100_seed79_latticeseed21_slowcool_spincorrelationfunctiontot
# fcc8x8x4yopen_nnJyz1p04_nnnJy0p67_T30K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot

#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x8yopen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_T1p5K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Jyz1_Jy0p5/fcc8x8x8yopen_nnJyz1p04_nnJxy0p3_nnnJy0p67_T15K_eq10000_mc10000_bins100_seed79_latticeseed21_II_slowcool_spincorrelationfunctiontot.txt", "r")

# The October files:
# DON'T have nnnJy0p67 in the files in the line below. That just remained in the name by accident.
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc8x8x8yopen_T5_nnJyz1p04_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_seed59_latticeseed23_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc8x8x8yopen_nnJxy0p2977_nnJyz1p036_nnnJy0p6701_T30K_eqst100000_mc10000_bins100_latticeseed31_MCseed71_slowcool_date221017_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc8x8x8yopen_nnJxy0p2977_nnJyz1p036_nnJxzm0p1121_nnnJy0p6701_nnnJz0p0469_sianDx0p3392_Dy1p8194_T25K_eqst100000_mc10000_bins100_latticeseed31_MCseed71_slowcool_date221017_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc14x14x14yopen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p7_nnnJzm0p05_sianDx0p34_Dy1p82_T10K_eq100000_mc100000_bins100_seed31_latticeseed71_slowcool_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc10x10x10yopen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p7_nnnJzm0p05_sianDx0p34_Dy1p82_T30K_eq100000_mc10000_bins100_seed31_latticeseed71_slowcool_spincorrelationfunctiontot.txt", "r")

# These were named wrong. It says 14x14x14, but it is really 20x20x20. Plus, T=10K has a slightly different name for both
# 30.10.17:
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc14x14x14yopen_nnJxy0p01_nnJyz1p04_nnnJy0p6701_incommtest_10K_eqst100000_mc10000_bins100_latticeseed31_MCseed71_slowcool_date301017_spincorrelationfunctiontot.txt", "r")
# 31.10.17:
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc14x14x14yopen_nnJxy0p1_nnJyz1p04_nnnJy0p6701_incommtest_T30K_eqst100000_mc10000_bins100_latticeseed31_MCseed71_slowcool_date311017_spincorrelationfunctiontot.txt", "r")
#infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc20x20x20yopen_nnJyz1p3337_nnnJy0p6701_incommtest_T30K_eqst100000_mc10000_bins100_latticeseed31_MCseed71_slowcool_date311017_spincorrelationfunctiontot.txt", "r")
# 03.11.17:
infile = open("/home/ubu/Master/Spincorrelation_Results/Frustration_fcc/Data_Jyz1p04_Jy0p67/fcc20x20x20yopen_nnJyz1p3337_nnnJy0p6701_incommtest_T12p5K_eqst100000_mc10000_bins100_latticeseed31_MCseed71_slowcool_date031117_spincorrelationfunctiontot.txt", "r")

# FOR HEAVEN'S SAKE, REMEMBER TO UPDATE THIS:

L1 = 20
L2 = 20
L3 = 20

# Getting the indices of the sites along (0,y,0)
yinfile = open("/home/ubu/Master/Spincorrelation_Results/Yline_0y0/fcc20x20x20_qyline.txt", "r")
# Getting lists ready to store the data
ns  = []

# Read the lines
ylines = yinfile.readlines()

# Getting data from the file
for yline in ylines:
    words = yline.split()
    if len(words) != 0:
        n = int(words[0])
        ns.append(n)
        
# Remember to close the file
yinfile.close()


firstline = infile.readline() # Reads the first line only # Should maybe expand this. But who knows
beta, N = firstline.split()

beta = float(beta); N = float(N)

T = 11.6045221/beta

# Getting lists ready to store the data
spin        = []
spinstd     = []

# Read the lines
lines = infile.readlines()

# Getting data from the file
for line in lines:
    words = line.split()
    if len(words) != 0:
        spinput = float(words[0])
        spin.append(spinput) 
        spinput = float(words[1])
        spinstd.append(spinput)
        
# Remember to close the file
infile.close()

# Converting to arrays
spin    = array(spin)
spinstd = array(spinstd)
ns      = array(ns)
Kfactor = 2./L1 # With our implementation so far. Should do something smarter to allow L1!=L2!=L3. # Look it up
# I should rather extract something. Well, this will do for now

length = len(ns)

ylinespins     = zeros(length) 
ylinespins_std = zeros(length)
Ks             = zeros(length)

for i in range(0,length):    # There should be a better way to do this. I guess we want to plot against beta...
    index = ns[i]
    ylinespins[i]     = spin[index]
    ylinespins_std[i] = spinstd[index]
    Ks[i] = Kfactor*i
    print "Spin correlation function: ", ylinespins[i], " standard deviation: ", ylinespins_std[i]

print ylinespins[length/2]

### Plotting
# If I want errorbars

figure(figsize=(6,5))
#ticklabel_format(style='sci', axis='y', scilimits=(0,0))
errorbar(Ks, ylinespins, yerr=ylinespins_std, capsize=2)  #, fmt=None) #Will maybe want to turn fmt off again...
title(r'Spin correlation function for T=%.2f K, %ix%ix%i fcc'% (T, L1, L2, L3) , fontsize=16)
xlabel(r'$K$', fontsize=20)
ylabel(r'$<\vec{S}_\vec{q}\vec{S}_{-\vec{q}}>$', fontsize=20)
axis([-0.1, max(Ks)*1.01, (min(ylinespins)-max(ylinespins_std))*0.9, (max(ylinespins)+max(ylinespins_std))*1.01])
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
show()


#maxspin = max(ylinespins)


'''
figure()
plot(qxs, spin, 'o')  #, fmt=None) #Will maybe want to turn fmt off again...
title('Spin correlation function for beta=%.2f, %i p fcc'% (beta, N) , fontsize=16)
xlabel(r'$q_x$', fontsize=20)
ylabel(r'$S^z_q$', fontsize=20)
axis([-0.0000001, maxqx, -0.0000001, 1.01*maxspin])
show()
'''
     
  
