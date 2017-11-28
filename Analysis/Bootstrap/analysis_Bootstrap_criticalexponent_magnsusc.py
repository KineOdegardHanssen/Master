from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys
from random import randint
from scipy.stats import norm
from scipy.optimize import curve_fit

'''
# Test of rng
for i in range(8):
    print randint(0,3)
'''

# I will extract only info about Tc from the Binder cumulant first

def giveline(a,b,x):
    return a*x+b

def givequadratic(a, b, c, x): # Should this really be quadratic? I look at a small interval...
    return a*x**2+b*x+c
    
def givecubic(a,b,c,d,x):
    return a*x**3+b*x**2+c*x+d

def givequadruple(a,b,c,d,e,x):
    return a*x**4+b*x**3+c*x**2+d*x+e
    
def givepowerlaw(k,n,x):
    return k*x**n

def gaus(x,a,x0,sigma, d):
    return a*exp(-(x-x0)**2/(2*sigma**2))+d

def findmaxima(beta, f1, warning):
    arraylength = len(f1)
    maxholder = -10000  # Definitely smaller than our maxima
    betaatmax = -10000  # So we will react when something goes wrong. In case the warning isn't enough.
    indexofmax = -10000 # Same
    for i in range(arraylength):
        if f1[i]>maxholder: # Just a safeguard. We're doomed if this happens anyway
            maxholder = f1[i]
            betaatmax = beta[i]
            indexofmax = i
    if warning<3:
        if indexofmax==(arraylength-1):
            print "WARNING! The largest element is also the last!"
    return betaatmax, maxholder

def fittingfunc(x, crossing, weight, exponent): # x = L**(-1)
    return crossing + weight*x**(exponent)

def give_fittingfunc(x, crossing, weight, exponent): # x = L**(-1)
    lenx = len(x)
    expfunc = zeros(lenx)
    for i in range(lenx):
        expfunc[i] = crossing + weight*x[i]**(exponent)
    return expfunc

def converttemp(crossingbeta):
    return 11.6045221/crossingbeta

def bootstrapping_each(no_of_betasA, binofbeta_firstindexA, binofbeta_lastindexA, nbinsA, mxsc_bavsA, mysc_bavsA, mzsc_bavsA, mxsc_abs_bavsA, mysc_abs_bavsA, mzsc_abs_bavsA, mxssqc_bavsA, myssqc_bavsA, mzssqc_bavsA, mxsquadc_bavsA, mysquadc_bavsA, mzsquadc_bavsA):
    magnsuscsx = []; magnsuscsy = []; magnsuscsz = []
    heatcaps = []
    for i in range(0, no_of_betasA): # For doing it for every temperature           
        nbinsthis = int(nbinsA[i])
        ### Reset all quantities I want to find by bootstrap    
        mx = 0; mxabs = 0; mx2 = 0; mx4 = 0; my = 0; myabs = 0; my2 = 0; my4 = 0; mz = 0; mzabs = 0; mz2 = 0; mz4 = 0;
        eav = 0; esqav = 0; hcav = 0;
        for j in range(0, int(nbinsA[i])): # For finding the averages
            ###Draw a random integer in [binofbeta_firstindexA[i], binofbeta_lastindexA[i]].
            n = randint(binofbeta_firstindexA[i], binofbeta_lastindexA[i])  # Draw a random number n times
            ###Extract O(n), add to average function
            eav += energy_bavsA[n]; esqav += energy_sq_bavsA[n]; hcav += hc_bavsA[n]
            mx += mxsc_bavsA[n]; mxabs += mxsc_abs_bavsA[n]; mx2 += mxssqc_bavsA[n]; mx4 += mxsquadc_bavsA[n];
            my += mysc_bavsA[n]; myabs += mysc_abs_bavsA[n]; my2 += myssqc_bavsA[n]; my4 += mysquadc_bavsA[n];
            mz += mzsc_bavsA[n]; mzabs += mzsc_abs_bavsA[n]; mz2 += mzssqc_bavsA[n]; mz4 += mzsquadc_bavsA[n];
        ### Divide the average of O(n) by no_of_bins_each_beta to ACTUALLY make the average
        eav = eav/nbins; esqav = esqav/nbins; hcav = hcav/nbins
        mx = mx/nbins; mxabs = mxabs/nbins; mx2 = mx2/nbins; mx4 = mx4/nbins;
        my = my/nbins; myabs = myabs/nbins; my2 = my2/nbins; my4 = my4/nbins;
        mz = mz/nbins; mzabs = mzabs/nbins; mz2 = mz2/nbins; mz4 = mz4/nbins;
        ### Probably store it in a list. Ya know, for standard deviations and such
        ### Magnetic susceptibility per spin
        magnsuscx_this = betasA[i]*(mx2-(mxabs*mxabs))*N # OFS said I should use the absolute magn to the first power. AF, prob 
        magnsuscy_this = betasA[i]*(my2-(myabs*myabs))*N
        magnsuscz_this = betasA[i]*(mz2-(mzabs*mzabs))*N
        ### Total heat capacity
        heatcap_this = betasA[i]**2*(esqav-(eav*eav))/N
        ### Feed into lists
        magnsuscsx.append(magnsuscx_this)
        magnsuscsy.append(magnsuscy_this)
        magnsuscsz.append(magnsuscz_this)
        heatcaps.append(heatcap_this)
    return magnsuscsx, magnsuscsy, magnsuscsz, heatcaps 


def max_temp(LA, filenameA, quadraticnotcubicfit, toplotornot, usegaussian, cutit, cutlower, cutupper): # Extract the beta value of the crossing between two graphs
    N = LA*LA*LA
    infileA = open(filenameA, "r")
    # Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
    
    # Sorting arrays
    betasA                 = []          # List of beta values
    binofbeta_firstindexA  = []          # I guess nested lists are a bit of a hassle
    binofbeta_lastindexA   = []          # Could also just feed the number of bins in. But, the less such input, the better, probably
    # Is the above unsafe? Should I just find the no of bins and add instead, a la, i*nbins+something?
    
    # Other quantities
    energy_bavsA      = []          # Index 1
    energy_sq_bavsA   = []          # Index 2
    hc_bavsA          = []          # Index 3
    mxsc_bavsA        = []          # Index 4
    mxsc_abs_bavsA    = []          # Index 5
    mxssqc_bavsA      = []          # Index 6
    mxsquadc_bavsA    = []          # Index 7
    mysc_bavsA        = []          # Index 8
    mysc_abs_bavsA    = []          # Index 9
    myssqc_bavsA      = []          # Index 10
    mysquadc_bavsA    = []          # Index 11
    mzsc_bavsA        = []          # Index 12
    mzsc_abs_bavsA    = []          # Index 13
    mzssqc_bavsA      = []          # Index 14
    mzsquadc_bavsA    = []          # Index 15
    
    # Don't think I need the magnetizations of q=(0,0,0), but we'll see...
    # Read the rest of the lines
    lines = infileA.readlines()  # This works, I have checked it
    
    i = 0 # Some sort of counter to get the indices right
    betabefore = 0 # Test the betas so we can gather the bins.
    binofbeta_firstindexA.append(0) 
    # Should I have some list over the indices of the betas belonging to each bin?
    
    # Getting data from the file
    for line in lines:
        words = line.split()
        if len(words) != 0:
            # Betas
            beta = float(words[0])
            if beta != betabefore:
                betasA.append(beta)
                betabefore = beta
                if i != 0:
                    binofbeta_firstindexA.append(i)
                    binofbeta_lastindexA.append(i-1)
            # energy
            en = float(words[1])
            energy_bavsA.append(en)          # Index 1
            # energy squared
            en = float(words[2])
            energy_sq_bavsA.append(en)       # Index 2
            # heat capacity
            en = float(words[3])
            hc_bavsA.append(en)              # Index 3
            # <m_x(q)>
            en = float(words[4])
            mxsc_bavsA.append(en)            # Index 4
            # <|m_x(q)|>
            en = float(words[5])
            mxsc_abs_bavsA.append(en)        # Index 5
            # <m^2_x(q)>
            en = float(words[6])
            mxssqc_bavsA.append(en)          # Index 6
            # <m^4_x(q)>
            en = float(words[7])
            mxsquadc_bavsA.append(en)        # Index 7
            # <m_y(q)>
            en = float(words[8])
            mysc_bavsA.append(en)            # Index 8
            # <|m_y(q)|>
            en = float(words[9])
            mysc_abs_bavsA.append(en)        # Index 9
            # <m^2_y(q)>
            en = float(words[10])
            myssqc_bavsA.append(en)          # Index 10
            # <m^4_y(q)>
            en = float(words[11])
            mysquadc_bavsA.append(en)        # Index 11
            # <m_z(q)>
            en = float(words[12])
            mzsc_bavsA.append(en)            # Index 12
            # <|m_z(q)|>
            en = float(words[13])
            mzsc_abs_bavsA.append(en)        # Index 13
            # <m^2_z(q)>
            en = float(words[14])
            mzssqc_bavsA.append(en)          # Index 14
            # <m^4_z(q)>
            en = float(words[15])
            mzsquadc_bavsA.append(en)        # Index 15
        i += 1 # Increasing the counter
        
    binofbeta_lastindexA.append(i-1)  # The index has been updated one time too many
    
    # Remember to close the file
    infileA.close()
    
    nbinsA = zeros(len(betasA))
    for i in range(0,len(betasA)):
        nbinsA[i] = binofbeta_lastindexA[i]-binofbeta_firstindexA[i]+1
    
    no_of_bins_each_betaA = binofbeta_lastindexA[0]+1 # Due to the indexing in Python + we don't allow
    no_of_betasA          = len(betasA)               # The number of different temperatures ### Should I convert to T now?
    no_of_BootstraprunsA  = 300                       # Good enough?
    
    betasA = array(betasA) # Not sure if this is neccessary, but maybe...
    
    if cutit==1:   # If we have simulations for a larger temperature interval than we wish
        betamin = cutlower
        betamax = cutupper
        # For A
        smallbetaarray = []
        binofbeta_firstindex_thosewewant = []
        binofbeta_lastindex_thosewewant = []
        nbins_wewant = []
        no_of_betas_wewant= 0
        for i in range(0, no_of_betasA):
            if (betasA[i]<betamax and betasA[i]>betamin):
                smallbetaarray.append(betasA[i])
                nbins_wewant.append(nbinsA[i])
                binofbeta_firstindex_thosewewant.append(binofbeta_firstindexA[i])
                binofbeta_lastindex_thosewewant.append(binofbeta_lastindexA[i])
                no_of_betas_wewant+=1
        no_of_betasA = no_of_betas_wewant
        betasA = smallbetaarray
        nbinsA = nbins_wewant
        binofbeta_firstindexA = binofbeta_firstindex_thosewewant
        binofbeta_lastindexA = binofbeta_lastindex_thosewewant
    betasA = array(betasA) # Not sure if this is neccessary, but maybe...
    
    if toplotornot==0:
        plotter = 1
    else:
        plotter = 0
    warning = 0
    betaatmax_av = 0
    betas_at_max = zeros(no_of_BootstraprunsA)
    maxtemp_av = 0
    maxtemps = zeros(no_of_BootstraprunsA)
    tempatmax_K_av = 0
    temp_at_max = zeros(no_of_BootstraprunsA)
    for k in range(0, no_of_BootstraprunsA): # Want to repeat the random selection of bins a number of times
        # Okay, maybe I don't need to make a call  
        #magnsuscsx, magnsuscsy, magnsuscsz, heatcaps = bootstrapping_each(no_of_betasA, binofbeta_firstindexA, binofbeta_lastindexA, nbinsA, mxsc_bavsA, mysc_bavsA, mzsc_bavsA, mxsc_abs_bavsA, mysc_abs_bavsA, mzsc_abs_bavsA, mxssqc_bavsA, myssqc_bavsA, mzssqc_bavsA, mxsquadc_bavsA, mysquadc_bavsA, mzsquadc_bavsA)
        magnsuscsx = []; magnsuscsy = []; magnsuscsz = []
        heatcaps = []
        for i in range(0, no_of_betasA): # For doing it for every temperature           
            nbins = int(nbinsA[i])
            ### Reset all quantities I want to find by bootstrap    
            mx = 0; mxabs = 0; mx2 = 0; mx4 = 0; my = 0; myabs = 0; my2 = 0; my4 = 0; mz = 0; mzabs = 0; mz2 = 0; mz4 = 0;
            eav = 0; esqav = 0; hcav = 0;
            for j in range(0, nbins): # For finding the averages
                ###Draw a random integer in [binofbeta_firstindexA[i], binofbeta_lastindexA[i]].
                n = randint(binofbeta_firstindexA[i], binofbeta_lastindexA[i])  # Draw a random number n times
                ###Extract O(n), add to average function
                eav += energy_bavsA[n]; esqav += energy_sq_bavsA[n]; hcav += hc_bavsA[n]
                mx += mxsc_bavsA[n]; mxabs += mxsc_abs_bavsA[n]; mx2 += mxssqc_bavsA[n]; mx4 += mxsquadc_bavsA[n];
                my += mysc_bavsA[n]; myabs += mysc_abs_bavsA[n]; my2 += myssqc_bavsA[n]; my4 += mysquadc_bavsA[n];
                mz += mzsc_bavsA[n]; mzabs += mzsc_abs_bavsA[n]; mz2 += mzssqc_bavsA[n]; mz4 += mzsquadc_bavsA[n];
            ### Divide the average of O(n) by no_of_bins_each_beta to ACTUALLY make the average
            eav = eav/nbins; esqav = esqav/nbins; hcav = hcav/nbins
            mx = mx/nbins; mxabs = mxabs/nbins; mx2 = mx2/nbins; mx4 = mx4/nbins;
            my = my/nbins; myabs = myabs/nbins; my2 = my2/nbins; my4 = my4/nbins;
            mz = mz/nbins; mzabs = mzabs/nbins; mz2 = mz2/nbins; mz4 = mz4/nbins;
            ### Probably store it in a list. Ya know, for standard deviations and such
            ### Magnetic susceptibility per spin
            magnsuscx_this = betasA[i]*(mx2-(mxabs*mxabs))*N # OFS said I should use the absolute magn to the first power. AF, prob 
            magnsuscy_this = betasA[i]*(my2-(myabs*myabs))*N
            magnsuscz_this = betasA[i]*(mz2-(mzabs*mzabs))*N
            ### Total heat capacity
            heatcap_this = betasA[i]**2*(esqav-(eav*eav))/N
            ### Feed into lists
            magnsuscsx.append(magnsuscx_this)
            magnsuscsy.append(magnsuscy_this)
            magnsuscsz.append(magnsuscz_this)
            heatcaps.append(heatcap_this)       
        ### Do a fitting of the results. To a line, a quadr. func., a pol. of the third degree. Check this. OFS suggested quadr
        if quadraticnotcubicfit==0:
            magnsuscsz = array(magnsuscsz) 
            fitvectorzcompA = polyfit(betasA, magnsuscsz, 1) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1];
                        
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = giveline(azA, bzA, fbetas) # Should this really be quadratic? I look at a small interval...
        if quadraticnotcubicfit==1:
            magnsuscsz = array(magnsuscsz) 
            fitvectorzcompA = polyfit(betasA, magnsuscsz, 2) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givequadratic(azA, bzA, czA, fbetas)
        if quadraticnotcubicfit==2:
            # For A
            magnsuscsz = array(magnsuscsz) 
            fitvectorzcompA = polyfit(betasA, magnsuscsz, 3) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]; dzA = fitvectorzcompA[3];
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givecubic(azA, bzA, czA, dzA, fbetas)
        if quadraticnotcubicfit==3:
            # For A
            magnsuscsz = array(magnsuscsz) 
            fitvectorzcompA = polyfit(betasA, magnsuscsz, 4) # Fits the function points to a quadruple polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]; dzA = fitvectorzcompA[3];ezA = fitvectorzcompA[4];
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givequadruple(azA, bzA, czA, dzA, ezA, fbetas) # Should this really be quadruple? I look at a small interval...
        if usegaussian!=0:
            n = len(betasA)
            #mean = sum(magnsuscsz*betasA)/n
            mean, thing = findmaxima(betasA, magnsuscsz, warning)
            sigma = sum(magnsuscsz*(betasA-mean)**2)/n
            d = min(magnsuscsz)
            popt,pcov = curve_fit(gaus,betasA,magnsuscsz,p0=[1,mean,sigma,d])
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            fA = gaus(fbetas,*popt)
            #print "In Gaussian fit!"
        #print "magnsuscsz[0] =", magnsuscsz[0],"magnsuscsz[1] =", magnsuscsz[1]

        betaatmax, maxtemp = findmaxima(fbetas, fA, warning)
        betaatmax_av += betaatmax
        betas_at_max[k] = betaatmax
        maxtemp_av += maxtemp
        maxtemps[k] = maxtemp
        tmK = converttemp(betaatmax)
        tempatmax_K_av += tmK
        temp_at_max[k] = tmK
        if warning<3:
            warning += 1
        
        if plotter==0:
            print "beta at max =", betaatmax, "; temp at max =", tmK
            figure(figsize=(6,5))
            plot(fbetas, fA, label='Fit')
            hold('on')
            plot(betasA, magnsuscsz, 'ro', label='Output')
            title(r'Magnetic susceptibility for L=%i'%LA, fontsize=14)
            xlabel(r'$\beta$', fontsize=20)
            ylabel(r'$\chi_{L,z}$', fontsize=20)
            legend(loc="lower center")
            tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
            show()
            if k==2: # Getting more plots
                plotter = 1
        
    betaatmax_av = betaatmax_av/no_of_BootstraprunsA
    betaatmax_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        betaatmax_stddev += (betaatmax_av-betas_at_max[i])*(betaatmax_av-betas_at_max[i])
    betaatmax_stddev = sqrt(betaatmax_stddev/(no_of_BootstraprunsA-1))
    
    maxtemp_av = maxtemp_av/no_of_BootstraprunsA
    maxtemp_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        maxtemp_stddev += (maxtemp_av-maxtemps[i])*(maxtemp_av-maxtemps[i])
    maxtemp_stddev = sqrt(maxtemp_stddev/(no_of_BootstraprunsA-1))

    tempatmax_K_av = tempatmax_K_av/no_of_BootstraprunsA
    tempatmax_K_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        tempatmax_K_stddev += (tempatmax_K_av-temp_at_max[i])*(tempatmax_K_av-temp_at_max[i])
    tempatmax_K_stddev = sqrt(tempatmax_K_stddev/(no_of_BootstraprunsA-1))
    
    return betaatmax_av, betaatmax_stddev, maxtemp_av, maxtemp_stddev, tempatmax_K_av, tempatmax_K_stddev

#################### This is where we give our input #######################
LA = 6
LB = 8
LC = 10
LD = 12
LE = 14
LF = 16

#### 1 000 000 MCsteps/bin
# Jensen
'''
filenameA = "fcc6x6x6yopen_verylongrunformagnsusc_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameB = "fcc8x8x8yopen_verylongrunformagnsusc_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameC = "fcc10x10x10yopen_verylongrunformagnsusc_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameD = "fcc12x12x12yopen_verylongrunformagnsusc_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameE = "fcc14x14x14yopen_verylongrunformagnsusc_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameF = "fcc16x16x16yopen_verylongrunformagnsusc_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
'''
#### 100 000 MCsteps/bin
# Jensen
filenameA = "fcc6x6x6yopen_magnsusc_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_rightinterval_binavgs_better.txt"
filenameB = "fcc8x8x8yopen_magnsusc_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_rightinterval_binavgs_better.txt"
filenameC = "fcc10x10x10yopen_magnsusc_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_rightinterval_binavgs_better.txt"
filenameD = "fcc12x12x12yopen_magnsusc_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_rightinterval_binavgs_better.txt"
filenameE = "fcc14x14x14yopen_magnsusc_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_rightinterval_binavgs_better.txt"
filenameF = "fcc16x16x16yopen_magnsusc_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_rightinterval_binavgs_better.txt"

##################                                   #######################
bool16 = 1 #0 - we have only up to L=14 / anythong else- we have L=16 too

quadraticnotcubicfit = 3 # 0 gives a linear fit to the Binder cumulants, 1 gives a quadratic fit, 2 gives cubic, 3 gives quadruple
toplotornot = 1 # 0, don't plot, anything else leads to a plot
usegaussian = 0 # 0, don't, anything else leads to a Gaussian fitting

# If we only want to look at a small temperature interval
cutit = 0 # 0: Find a fit for the whole beta range; 1: Find a fit for a small interval

cutlower6 = 0.806
cutupper6 = 0.821

cutlower8 = 0.79
cutupper8 = 0.815

cutlower10 = 0.7
cutupper10 = 0.815

cutlower12 = 0.7
cutupper12 = 0.806

cutlower14 = 0.777
cutupper14 = 0.801

cutlower16 = 0.782
cutupper16 = 0.798

cutlower1 = 0.8075
#cutupper1 = 0.84
cutupper1 = 0.8225

bm_avA, bm_stddevA, tmax_avA, tmax_stddevA, tempmax_avA, tempmax_stddevA = max_temp(LA, filenameA, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower6, cutupper6)
bm_avB, bm_stddevB, tmax_avB, tmax_stddevB, tempmax_avB, tempmax_stddevB = max_temp(LB, filenameB, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower8, cutupper8)
bm_avC, bm_stddevC, tmax_avC, tmax_stddevC, tempmax_avC, tempmax_stddevC = max_temp(LC, filenameC, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower10, cutupper10)
bm_avD, bm_stddevD, tmax_avD, tmax_stddevD, tempmax_avD, tempmax_stddevD = max_temp(LD, filenameD, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower12, cutupper12)
bm_avE, bm_stddevE, tmax_avE, tmax_stddevE, tempmax_avE, tempmax_stddevE = max_temp(LE, filenameE, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower14, cutupper14)
if bool16!=0:
    bm_avF, bm_stddevF, tmax_avF, tmax_stddevF, tempmax_avF, tempmax_stddevF = max_temp(LF, filenameF, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower14, cutupper14)

betaatmaxes      = [bm_avA, bm_avB, bm_avC, bm_avD, bm_avE]
betaatmaxes_stdv = [bm_stddevA, bm_stddevB, bm_stddevC, bm_stddevD, bm_stddevE]
tempatmaxes      = [tempmax_avA, tempmax_avB, tempmax_avC, tempmax_avD, tempmax_avE]
tempatmaxes_stdv = [tempmax_stddevA, tempmax_stddevB, tempmax_stddevC, tempmax_stddevD, tempmax_stddevE]
tmaxes           = [tmax_avA, tmax_avB, tmax_avC, tmax_avD, tmax_avE]
tmaxes_stdv      = [tmax_stddevA, tmax_stddevB, tmax_stddevC, tmax_stddevD, tmax_stddevE]
vec1dL           = [1./LA, 1./LB, 1./LC, 1./LD, 1./LE]

if bool16!=0:
    betaatmaxes.append(bm_avF)
    betaatmaxes_stdv.append(bm_stddevF)
    tempatmaxes.append(tempmax_avF)
    tempatmaxes_stdv.append(tempmax_stddevF)
    tmaxes.append(tmax_avF)
    tmaxes_stdv.append(tmax_stddevF)
    vec1dL.append(1./LF)

betaatmaxes      = array(betaatmaxes)
betaatmaxes_stdv = array(betaatmaxes_stdv)
tempatmaxes      = array(tempatmaxes)
tempatmaxes_stdv = array(tempatmaxes_stdv)
tmaxes           = array(tmaxes)
tmaxes_stdv      = array(tmaxes_stdv)
vec1dL           = array(vec1dL)

# Print to file 
outfile = open("maxmagnsusc_Jensen_006.txt","w")
for i in range(len(tmaxes)):
    outfile.write('%.16f %.16f \n'%(tmaxes[i], tmaxes_stdv[i]))
outfile.close()

outfileT = open("temp_at_maxmagnsusc_Jensen_006.txt","w")
for i in range(len(tmaxes)):
    outfileT.write('%.16f %.16f \n'%(tempatmaxes[i], tempatmaxes_stdv[i]))
outfileT.close()


max1dL_inplot = (max(vec1dL)*1.02)

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL, tmaxes, yerr=tmaxes_stdv, fmt=None, capsize=2)
hold('on')
title(r'Maxes of the magnetic susceptibility', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\chi_{max}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, max1dL_inplot, (min(tmaxes)-max(tmaxes_stdv))*0.99, (max(tmaxes)+max(tmaxes_stdv))*1.01])
show()

figure(figsize=(6,5))
errorbar(vec1dL, tempatmaxes, yerr=tempatmaxes_stdv, fmt=None, capsize=2)
hold('on')
title(r'Maxes of the magnetic susceptibility', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$T_{max}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, max1dL_inplot, (min(tempatmaxes)-max(tempatmaxes_stdv))*0.999, (max(tempatmaxes)+max(tempatmaxes_stdv))*1.001])
show()
