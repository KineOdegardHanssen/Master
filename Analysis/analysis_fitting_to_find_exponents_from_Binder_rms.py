from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys

def findline(Lbefl, Llast, ybefl, ylast): # Weird name to remember what goes where
    x0 = Lbefl; x1 = Llast; y0=ybefl; y1 = ylast;
    print "y0=%.3f, y1=%.3f"%(y0,y1)
    a = (y1-y0)/float(x1-x0) # Float just in case
    b = y1-(a*x1)
    print "a=%.3f, b=%.3f"%(a,b)
    return a, b

def makeline(x,a,b):
    return x*a+b

#Getting the input data
# L, L+2
infileT = open("LandLp2_crtemps_data_Li_001.txt","r")
infileU = open("LandLp2_crUs_data_Li_001.txt","r")
infileT4 = open("LandLp4_crtemps_data_Li_001.txt","r")
infileU4 = open("LandLp4_crUs_data_Li_001.txt","r")
lines = infileT.readlines()

crT      = []
crT_stdv = []

for line in lines:
    words = line.split()
    if len(words) != 0:
        crT.append(float(words[0]))
        crT_stdv.append(float(words[1])) 

infileT.close()

lines = infileU.readlines()

crU      = []
crU_stdv = []

for line in lines:
    words = line.split()
    if len(words) != 0:
        crU.append(float(words[0]))
        crU_stdv.append(float(words[1])) 

infileU.close()

crT = array(crT) 
crT_stdv = array(crT_stdv)
crU = array(crU)
crU_stdv = array(crU_stdv)

# These are always the same for us now.
vec1dL = array([1./6, 1./8, 1./10, 1./12, 1./14])

# L, L+4


lines = infileT4.readlines()

crT4      = []
crT_stdv4 = []

for line in lines:
    words = line.split()
    if len(words) != 0:
        crT4.append(float(words[0]))
        crT_stdv4.append(float(words[1])) 

infileT4.close()



lines = infileU4.readlines()

crU4      = []
crU_stdv4 = []

for line in lines:
    words = line.split()
    if len(words) != 0:
        crU4.append(float(words[0]))
        crU_stdv4.append(float(words[1])) 

infileU4.close()

crT4 = array(crT4) 
crT_stdv4 = array(crT_stdv4)
crU4 = array(crU4)
crU_stdv4 = array(crU_stdv4)

vec1dL4 = array([1./6, 1./8, 1./10, 1./12])

# Plotting parameters:
T_min = (min(crT)-max(crT_stdv))*0.9999
T_max = (max(crT)+max(crT_stdv))*1.0001
U_min = (min(crU)-max(crU_stdv))*0.9
U_max = (max(crU)+max(crU_stdv))*1.1

'''
# Plotting T
figure(figsize=(6,5))
errorbar(vec1dL, crT, yerr=crT_stdv, capsize=2)
hold('on')
xlabel(r'$L^{-1}$', fontsize=20)
ylabel(r'T[K]', fontsize=20)
title(r'Fit of T to f.s.s., L, L+2', fontsize=14)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([min(vec1dL)*0.9, max(vec1dL)*1.1, T_min, T_max])
show()
'''
### Max size of error bars:

print "Max error bar, L,L+2, T:", crT_stdv
print "Max error bar, L,L+4, T:", crT_stdv4
print "Max error bar, L,L+2, U:", crU_stdv
print "Max error bar, L,L+4, U:", crU_stdv4

aerrT2 = average(crT_stdv)
aerrT4 = average(crT_stdv4)
aerrU2 = average(crU_stdv)
aerrU4 = average(crU_stdv4)

print "Average error, L,L+2, T:", aerrT2
print "Average error, L,L+4, T:", aerrT4
print "Average error, L,L+2, U:", aerrU2
print "Average error, L,L+4, U:", aerrU4

print "Relative difference in the average error, T:", (aerrT2-aerrT4)/aerrT2
print "Relative difference in the average error, U:", (aerrU2-aerrU4)/aerrU2

### Doing the analysis/Guessing
ngT = 2 # The number of guesses
ngU = 2 # The number of guesses
exponentsT = linspace(3,7,ngT) # Omega+1/nu
exponentsU = linspace(3,5,ngU) # Omega

len2 = len(vec1dL)
len4 = len(vec1dL4)
la  = len2-1
la4 = len4-1
dfl = 1 # Distance from last point # This is important
maxy = 0
miny = 0
shift = -1 # Shift for L,L+2

midarrayL = array([ 1./8, 1./10, 1./12])
crTm = zeros(3)
crTm_stdv = zeros(3)
crUm = zeros(3) 
crUm_stdv = zeros(3)
for i in range(3):
    crTm[i] = crT[i+1]
    crTm_stdv[i] = crT_stdv[i+1]
    crUm[i] = crU[i+1]
    crUm_stdv[i] = crU_stdv[i+1]

# For L, L+2:

print "L,L+2:"

sp = 0 # skip this number of points
npoints = len(crT)-sp

for i in range(ngT):
    #print " "
    plotvecLT = zeros(len2)
    plotvecLTalt = zeros(len2-sp)
    crTshort = zeros(len2-sp)
    crT_stdv_short = zeros(len2-sp)
    eT = exponentsT[i]
    for j in range(len2):
        plotvecLT[j] = vec1dL[j]**eT
    for j in range(len2-sp):
        plotvecLTalt[j] = vec1dL[j+sp]**eT
        crTshort[j] = crT[j+sp]
        crT_stdv_short[j] = crT_stdv[j+sp]
    lcrTshort = max(crTshort)+max(crT_stdv_short)
    scrTshort = min(crTshort)-max(crT_stdv_short)
    maxy= lcrTshort*1.0001
    miny= scrTshort*0.9999
    fitvectorzcompA = polyfit(plotvecLTalt, crTshort, 1) # Fits the function points to a quadratic polynomial
    a = fitvectorzcompA[0]; b = fitvectorzcompA[1];
    straightline = makeline(plotvecLTalt,a, b) # Should this really be quadratic? I look at a small interval...
    # Only looking at the last crossings
    figure(figsize=(7,5))
    errorbar(plotvecLTalt, crTshort, yerr=crT_stdv_short, capsize=2)
    hold('on')
    plot(plotvecLTalt, straightline, '--', label='Straight line')
    xlabel(r'$L^{-%.3f}$'%eT, fontsize=20)
    ylabel(r'T[K]', fontsize=20)
    title(r'Fit of T to f.s.s., L, L+2', fontsize=14)
    tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axis([min(plotvecLTalt)*0.9, max(plotvecLTalt)*1.1, miny, maxy])
    legend(loc='upper left')
    show()
    # Plotting the rms of the innermost points as well
    plotvecLTalt = zeros(len2-sp-1)
    crTshort = zeros(len2-sp-1)
    crT_stdv_short = zeros(len2-sp-1)
    for j in range(len2-sp-1):
        plotvecLTalt[j] = vec1dL[j+sp+1]**eT
        crTshort[j] = crT[j+sp+1]
        crT_stdv_short[j] = crT_stdv[j+sp+1]
    lcrTshort = max(crTshort)+max(crT_stdv_short)
    scrTshort = min(crTshort)-max(crT_stdv_short)
    maxy= lcrTshort*1.0001
    miny= scrTshort*0.9999
    fitvectorzcompA = polyfit(plotvecLTalt, crTshort, 1) # Fits the function points to a quadratic polynomial
    a = fitvectorzcompA[0]; b = fitvectorzcompA[1];
    straightline = makeline(plotvecLTalt,a, b) # Should this really be quadratic? I look at a small interval...
    # Only looking at the last crossings
    figure(figsize=(7,5))
    errorbar(plotvecLTalt, crTshort, yerr=crT_stdv_short, capsize=2)
    hold('on')
    plot(plotvecLTalt, straightline, '--', label='Straight line')
    xlabel(r'$L^{-%.3f}$'%eT, fontsize=20)
    ylabel(r'T[K]', fontsize=20)
    title(r'Fit of T to f.s.s., L, L+2', fontsize=14)
    tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axis([min(plotvecLTalt)*0.9, max(plotvecLTalt)*1.1, miny, maxy])
    legend(loc='upper left')
    show()

for i in range(ngU):
    #print " "
    plotvecLU = zeros(len2)
    plotvecLUalt = zeros(len2-sp)
    crUshort = zeros(len2-sp)
    crU_stdv_short = zeros(len2-sp)
    eU = exponentsU[i]
    for j in range(len2):
        plotvecLU[j] = vec1dL[j]**eU
    for j in range(len2-sp):
        plotvecLUalt[j] = vec1dL[j+sp]**eU
        crUshort[j] = crU[j+sp]
        crU_stdv_short[j] = crU_stdv[j+sp]
    lcrUshort = max(crUshort)+max(crU_stdv_short)
    scrUshort = min(crUshort)-max(crU_stdv_short)
    maxy= lcrUshort*1.0001
    miny= scrUshort*0.9999
    fitvectorzcompA = polyfit(plotvecLUalt, crUshort, 1) # Fits the function points to a quadratic polynomial
    a = fitvectorzcompA[0]; b = fitvectorzcompA[1];
    straightline = makeline(plotvecLUalt,a, b) # Should this really be quadratic? I look at a small interval...
    figure(figsize=(7,5))
    errorbar(plotvecLUalt, crUshort, yerr=crU_stdv_short, capsize=2)
    hold('on')
    plot(plotvecLUalt, straightline, '--', label='Straight line')
    xlabel(r'$L^{-%.3f}$'%eU, fontsize=20)
    ylabel(r'U', fontsize=20)
    title(r'Fit of U to f.s.s., L, L+4, %i points'%npoints, fontsize=14)
    tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axis([min(plotvecLUalt)*0.9, max(plotvecLUalt)*1.1, miny, maxy])
    legend(loc='upper right')
    show()
    # Repeating to look at the last points
    plotvecLUalt = zeros(len2-sp-1)
    crUshort = zeros(len2-sp-1)
    crU_stdv_short = zeros(len2-sp-1)
    for j in range(len2-sp-1):
        plotvecLUalt[j] = vec1dL[j+sp+1]**eU
        crUshort[j] = crU[j+sp+1]
        crU_stdv_short[j] = crU_stdv[j+sp+1]
    lcrUshort = max(crUshort)+max(crU_stdv_short)
    scrUshort = min(crUshort)-max(crU_stdv_short)
    maxy= lcrUshort*1.0001
    miny= scrUshort*0.9999
    fitvectorzcompA = polyfit(plotvecLUalt, crUshort, 1) # Fits the function points to a quadratic polynomial
    a = fitvectorzcompA[0]; b = fitvectorzcompA[1];
    straightline = makeline(plotvecLUalt,a, b) # Should this really be quadratic? I look at a small interval...
    figure(figsize=(7,5))
    errorbar(plotvecLUalt, crUshort, yerr=crU_stdv_short, capsize=2)
    hold('on')
    plot(plotvecLUalt, straightline, '--', label='Straight line')
    xlabel(r'$L^{-%.3f}$'%eU, fontsize=20)
    ylabel(r'U', fontsize=20)
    title(r'Fit of U to f.s.s., L, L+2, %i points'%npoints, fontsize=14)
    tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axis([min(plotvecLUalt)*0.9, max(plotvecLUalt)*1.1, miny, maxy])
    legend(loc='upper right')
    show()

# For L, L+4:
# Plotting parameters:
T_min = (min(crT4)-max(crT_stdv4))*0.9990
T_max = (max(crT4)+max(crT_stdv4))*1.0001
U_min = (min(crU4)-max(crU_stdv4))*0.99
U_max = (max(crU4)+max(crU_stdv4))*1.01


#ngT = 21 # The number of guesses
#ngU = 21 # The number of guesses
#exponentsT = linspace(1,21,ngT) # Omega+1/nu
#exponentsU = linspace(1,21,ngU) # Omega
dfl = 2

ngT = 3 # The number of guesses
ngU = 3 # The number of guesses
exponentsT = linspace(3,7,ngT) # Omega+1/nu
exponentsU = linspace(3,5,ngU) # Omega

sp4 = 0 # skip this number of points
npoints = len(crT4)-sp4
print "L,L+4:"
for i in range(ngT):
    #print " "
    plotvecLT = zeros(len4)
    plotvecLTalt = zeros(len4-sp4)
    crTshort = zeros(len4-sp4)
    crT_stdv_short = zeros(len4-sp4)
    eT = exponentsT[i]
    for j in range(len4):
        plotvecLT[j] = vec1dL4[j]**eT
    for j in range(len4-sp4):
        plotvecLTalt[j] = vec1dL4[j+sp4]**eT
        crTshort[j] = crT4[j+sp4]
        crT_stdv_short[j] = crT_stdv4[j+sp4]
    lcrTshort = max(crTshort)+max(crT_stdv_short)
    scrTshort = min(crTshort)-max(crT_stdv_short)
    maxy= lcrTshort*1.0001
    miny= scrTshort*0.9999
    fitvectorzcompA = polyfit(plotvecLTalt, crTshort, 1) # Fits the function points to a quadratic polynomial
    a = fitvectorzcompA[0]; b = fitvectorzcompA[1];
    straightline = makeline(plotvecLTalt,a, b) # Should this really be quadratic? I look at a small interval...
    # Only looking at the last crossings
    figure(figsize=(7,5))
    errorbar(plotvecLTalt, crTshort, yerr=crT_stdv_short, capsize=2)
    hold('on')
    plot(plotvecLTalt, straightline, '--', label='Straight line')
    xlabel(r'$L^{-%.3f}$'%eT, fontsize=20)
    ylabel(r'T[K]', fontsize=20)
    title(r'Fit of T to f.s.s., L, L+4', fontsize=14)
    tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axis([min(plotvecLTalt)*0.9, max(plotvecLTalt)*1.1, miny, maxy])
    legend(loc='upper left')
    show()
    # Plotting the rms of the innermost points as well
    plotvecLTalt = zeros(len4-sp4-1)
    crTshort = zeros(len4-sp4-1)
    crT_stdv_short = zeros(len4-sp4-1)
    for j in range(len4-sp4-1):
        plotvecLTalt[j] = vec1dL4[j+sp4+1]**eT
        crTshort[j] = crT4[j+sp4+1]
        crT_stdv_short[j] = crT_stdv4[j+sp4+1]
    lcrTshort = max(crTshort)+max(crT_stdv_short)
    scrTshort = min(crTshort)-max(crT_stdv_short)
    maxy= lcrTshort*1.0001
    miny= scrTshort*0.9999
    fitvectorzcompA = polyfit(plotvecLTalt, crTshort, 1) # Fits the function points to a quadratic polynomial
    a = fitvectorzcompA[0]; b = fitvectorzcompA[1];
    straightline = makeline(plotvecLTalt,a, b) # Should this really be quadratic? I look at a small interval...
    # Only looking at the last crossings
    figure(figsize=(7,5))
    errorbar(plotvecLTalt, crTshort, yerr=crT_stdv_short, capsize=2)
    hold('on')
    plot(plotvecLTalt, straightline, '--', label='Straight line')
    xlabel(r'$L^{-%.3f}$'%eT, fontsize=20)
    ylabel(r'T[K]', fontsize=20)
    title(r'Fit of T to f.s.s., L, L+4', fontsize=14)
    tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axis([min(plotvecLTalt)*0.9, max(plotvecLTalt)*1.1, miny, maxy])
    legend(loc='upper left')
    show()

for i in range(ngU):
    #print " "
    plotvecLU = zeros(len4)
    plotvecLUalt = zeros(len4-sp4)
    crUshort = zeros(len4-sp4)
    crU_stdv_short = zeros(len4-sp4)
    eU = exponentsU[i]
    for j in range(len4):
        plotvecLU[j] = vec1dL4[j]**eU
    for j in range(len4-sp4):
        plotvecLUalt[j] = vec1dL4[j+sp4]**eU
        crUshort[j] = crU4[j+sp4]
        crU_stdv_short[j] = crU_stdv4[j+sp4]
    lcrUshort = max(crUshort)+max(crU_stdv_short)
    scrUshort = min(crUshort)-max(crU_stdv_short)
    maxy= lcrUshort*1.0001
    miny= scrUshort*0.9999
    fitvectorzcompA = polyfit(plotvecLUalt, crUshort, 1) # Fits the function points to a quadratic polynomial
    a = fitvectorzcompA[0]; b = fitvectorzcompA[1];
    straightline = makeline(plotvecLUalt,a, b) # Should this really be quadratic? I look at a small interval...
    figure(figsize=(7,5))
    errorbar(plotvecLUalt, crUshort, yerr=crU_stdv_short, capsize=2)
    hold('on')
    plot(plotvecLUalt, straightline, '--', label='Straight line')
    xlabel(r'$L^{-%.3f}$'%eU, fontsize=20)
    ylabel(r'U', fontsize=20)
    title(r'Fit of U to f.s.s., L, L+4, %i points'%npoints, fontsize=14)
    tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axis([min(plotvecLUalt)*0.9, max(plotvecLUalt)*1.1, miny, maxy])
    legend(loc='upper right')
    show()
    # Repeating to look at the last points
    plotvecLUalt = zeros(len4-sp4-1)
    crUshort = zeros(len4-sp4-1)
    crU_stdv_short = zeros(len4-sp4-1)
    for j in range(len4-sp4-1):
        plotvecLUalt[j] = vec1dL4[j+sp4+1]**eU
        crUshort[j] = crU4[j+sp4+1]
        crU_stdv_short[j] = crU_stdv4[j+sp4+1]
    lcrUshort = max(crUshort)+max(crU_stdv_short)
    scrUshort = min(crUshort)-max(crU_stdv_short)
    maxy= lcrUshort*1.0001
    miny= scrUshort*0.9999
    fitvectorzcompA = polyfit(plotvecLUalt, crUshort, 1) # Fits the function points to a quadratic polynomial
    a = fitvectorzcompA[0]; b = fitvectorzcompA[1];
    straightline = makeline(plotvecLUalt,a, b) # Should this really be quadratic? I look at a small interval...
    figure(figsize=(7,5))
    errorbar(plotvecLUalt, crUshort, yerr=crU_stdv_short, capsize=2)
    hold('on')
    plot(plotvecLUalt, straightline, '--', label='Straight line')
    xlabel(r'$L^{-%.3f}$'%eU, fontsize=20)
    ylabel(r'U', fontsize=20)
    title(r'Fit of U to f.s.s., L, L+4, %i points'%npoints, fontsize=14)
    tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    axis([min(plotvecLUalt)*0.9, max(plotvecLUalt)*1.1, miny, maxy])
    legend(loc='upper right')
    show()


minexpT = 4
mdlexpT = 5
maxexpT = 6

minomega = 3
mdlomega = 4
maxomega = 5

Lm3 = vec1dL**3
Lm4 = vec1dL**4
Lm5 = vec1dL**5
Lm6 = vec1dL**6

# Binder cumulant
abU3 = polyfit(Lm3, crU,1)
aU3 = abU3[0]; bU3 = abU3[1];

abU4 = polyfit(Lm4, crU,1)
aU4 = abU4[0]; bU4 = abU4[1];

abU5 = polyfit(Lm5, crU,1)
aU5 = abU5[0]; bU5 = abU5[1];

print "Binder cumulant:"
print "exp-stdv:"
print "weight = ", aU3, "; Uc = ", bU3

print "exp:"
print "weight = ", aU4, "; Uc = ", bU4

print "exp+stdv:"
print "weight = ", aU5, "; Uc = ", bU5

# temperature
abT4 = polyfit(Lm4, crT,1)
aT4 = abT4[0]; bT4 = abT4[1];

abT5 = polyfit(Lm5, crT,1)
aT5 = abT5[0]; bT5 = abT5[1];

abT6 = polyfit(Lm3, crT,1)
aT6 = abT6[0]; bT6 = abT6[1];

print "Temperature:"
print "exp-stdv:"
print "weight = ", aT4, "; Tc = ", bT4

print "exp:"
print "weight = ", aT5, "; Tc = ", bT5

print "exp-stdv:"
print "weight = ", aT6, "; Tc = ", bT6

