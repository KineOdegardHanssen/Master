from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys
from random import randint
from scipy.stats import norm
from scipy.optimize import curve_fit


def findline(x0, x1, y0, y1): # Weird name to remember what goes where
    #print "y0=%.3f, y1=%.3f"%(y0,y1)
    a = (y1-y0)/float(x1-x0) # Float just in case
    b = y1-(a*x1)
    #print "a=%.3f, b=%.3f"%(a,b)
    return a, b

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
    maxholder   = -10000  # Definitely smaller than our maxima
    maxvalue    = -10000  # So we will react when something goes wrong. In case the warning isn't enough.
    indexofmax  = -10000 # Same
    for i in range(arraylength):
        if f1[i]>maxholder: # Just a safeguard. We're doomed if this happens anyway
            maxholder = f1[i]
            maxvalue= beta[i]
            indexofmax = i
    if warning<3:
        if indexofmax==(arraylength-1):
            print "WARNING! The largest element is also the last!"
    return maxvalue, maxholder

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
            #if i==0:
            #    print en
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
    maxheatcap_av = 0
    maxheatcaps = zeros(no_of_BootstraprunsA)
    betaformax_av = 0
    betaformaxes = zeros(no_of_BootstraprunsA)
    tempformax_K_av = 0
    tempformaxes_K = zeros(no_of_BootstraprunsA)
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
            #if k==0:
            #    print eav
            ### Probably store it in a list. Ya know, for standard deviations and such
            heatcap_this = betasA[i]**2*(esqav-(eav*eav))/N
            heatcaps.append(heatcap_this)
        ### Do a fitting of the results. To a line, a quadr. func., a pol. of the third degree. Check this. OFS suggested quadr
        if quadraticnotcubicfit==0:
            heatcaps = array(heatcaps) 
            fitvectorzcompA = polyfit(betasA, heatcaps, 1) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1];
                        
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = giveline(azA, bzA, fbetas) # Should this really be quadratic? I look at a small interval...
        if quadraticnotcubicfit==1:
            heatcaps = array(heatcaps) 
            fitvectorzcompA = polyfit(betasA, heatcaps, 2) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givequadratic(azA, bzA, czA, fbetas)
        if quadraticnotcubicfit==2:
            heatcaps = array(heatcaps) 
            fitvectorzcompA = polyfit(betasA, heatcaps, 3) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]; dzA = fitvectorzcompA[3];
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givecubic(azA, bzA, czA, dzA, fbetas)
        if quadraticnotcubicfit==3:
            heatcaps = array(heatcaps) 
            fitvectorzcompA = polyfit(betasA, heatcaps, 4) # Fits the function points to a quadruple polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]; dzA = fitvectorzcompA[3];ezA = fitvectorzcompA[4];
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givequadruple(azA, bzA, czA, dzA, ezA, fbetas) # Should this really be quadruple? I look at a small interval...
        if usegaussian!=0:
            n = len(betasA)
            #mean = sum(magnsuscsz*betasA)/n
            mean, thing = findmaxima(betasA, heatcaps, warning)
            sigma = sum(heatcaps*(betasA-mean)**2)/n
            d = min(heatcaps)
            popt,pcov = curve_fit(gaus,betasA,heatcaps,p0=[1,mean,sigma,d])
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            fA = gaus(fbetas,*popt)
            #print "In Gaussian fit!"
        #print "magnsuscsz[0] =", magnsuscsz[0],"magnsuscsz[1] =", magnsuscsz[1]

        #for elem in range(len(fbetas)):
        #    print "fbeta(",elem,")= ", fbetas[elem]
        # Finding maxima
        betaformax, maxheatcap = findmaxima(fbetas, fA, warning)
        #print "betaformax:", betaformax, ", maxheatcap: ", maxheatcap
        # Max val, heat cap
        maxheatcap_av += maxheatcap
        maxheatcaps[k] = maxheatcap
        # Beta at max
        betaformax_av += betaformax
        betaformaxes[k] = betaformax
        # Temp in Kelvin at max
        tmK = converttemp(betaformax)
        tempformax_K_av += tmK
        tempformaxes_K[k] = tmK
        if warning<3:
            warning += 1
        
        if plotter==0:
            print "beta at max =", betaformax, "; temp at max =", tmK
            figure(figsize=(6,5))
            plot(fbetas, fA, label='Fit')
            hold('on')
            plot(betasA, heatcaps, 'ro', label='Output')
            title(r'Specific heat for L=%i'%LA, fontsize=14)
            xlabel(r'$\beta$', fontsize=20)
            ylabel(r'$c_{V}$', fontsize=20)
            legend(loc="lower center")
            tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
            show()
            if k==2: # Getting more plots
                plotter = 1

    maxheatcap_av = maxheatcap_av/no_of_BootstraprunsA
    maxheatcap_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        maxheatcap_stddev += (maxheatcap_av-maxheatcaps[i])*(maxheatcap_av-maxheatcaps[i])
    maxheatcap_stddev = sqrt(maxheatcap_stddev/(no_of_BootstraprunsA-1))
        
    betaformax_av = betaformax_av/no_of_BootstraprunsA
    betaformax_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        betaformax_stddev += (betaformax_av-betaformaxes[i])*(betaformax_av-betaformaxes[i])
    betaformax_stddev = sqrt(betaformax_stddev/(no_of_BootstraprunsA-1))

    tempformax_K_av = tempformax_K_av/no_of_BootstraprunsA
    tempformax_K_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        tempformax_K_stddev += (tempformax_K_av-tempformaxes_K[i])*(tempformax_K_av-tempformaxes_K[i])
    tempformax_K_stddev = sqrt(tempformax_K_stddev/(no_of_BootstraprunsA-1))
    
    return maxheatcap_av, maxheatcap_stddev, betaformax_av, betaformax_stddev, tempformax_K_av, tempformax_K_stddev

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
filenameA = "fcc6x6x6yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_binavgs.txt"
filenameB = "fcc8x8x8yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_binavgs.txt"
filenameC = "fcc10x10x10yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_binavgs.txt"
filenameD = "fcc12x12x12yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_binavgs.txt"
filenameE = "fcc14x14x14yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_binavgs.txt"
filenameF = "fcc16x16x16yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_slowcool_binavgs.txt"
'''
#### 100 000 MCsteps/bin
# Jensen
filenameA = "fcc6x6x6yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_binavgs.txt"
filenameB = "fcc8x8x8yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_binavgs_ri_2.txt"
filenameC = "fcc10x10x10yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_binavgs_ri_2.txt"
filenameD = "fcc12x12x12yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_binavgs_ri_new_better.txt"
filenameE = "fcc14x14x14yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_binavgs_ri.txt"
filenameF = "fcc16x16x16yopen_heatcap_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq100000_mc100000_bins100_divseeds_binavgs_ri.txt"

##################                                   #######################
bool16 = 1 #0 - we have only up to L=14 / anythong else- we have L=16 too

quadraticnotcubicfit = 3 # 0 gives a linear fit to the Binder cumulants, 1 gives a quadratic fit, 2 gives cubic, 3 gives quadruple
toplotornot = 1 # 0, don't plot, anything else leads to a plot
usegaussian = 0 # 0, don't, anything else leads to a Gaussian fitting

# If we only want to look at a small temperature interval
cutit = 0 # 0: Find a fit for the whole beta range; 1: Find a fit for a small interval
cutlower8 = 0.805
cutupper8 = 0.840

cutlower10 = 0.79
cutupper10 = 0.83

cutlower12 = 0.785
cutupper12 = 0.816

cutlower14 = 0.789
cutupper14 = 0.811

cutlower16 = 0.781
cutupper16 = 0.896

cutlower1 = 0.8075
cutupper1 = 0.84
cutupper1 = 0.8225
# Order: maxheatcap_av, maxheatcap_stddev, betaformax_av, betaformax_stddev, tempformax_K_av, tempformax_K_stddev
bm_avA, bm_stddevA, tmax_avA, tmax_stddevA, tempmax_avA, tempmax_stddevA = max_temp(LA, filenameA, quadraticnotcubicfit, toplotornot, usegaussian, cutit, cutlower1, cutupper1)
print "Done with A"
bm_avB, bm_stddevB, tmax_avB, tmax_stddevB, tempmax_avB, tempmax_stddevB = max_temp(LB, filenameB, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower8, cutupper8)
print "Done with B"
bm_avC, bm_stddevC, tmax_avC, tmax_stddevC, tempmax_avC, tempmax_stddevC = max_temp(LC, filenameC, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower10, cutupper10)
print "Done with C"
bm_avD, bm_stddevD, tmax_avD, tmax_stddevD, tempmax_avD, tempmax_stddevD = max_temp(LD, filenameD, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower12, cutupper12)
print "Done with D"
bm_avE, bm_stddevE, tmax_avE, tmax_stddevE, tempmax_avE, tempmax_stddevE = max_temp(LE, filenameE, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower14, cutupper14)
print "Done with E"
if bool16!=0:
    bm_avF, bm_stddevF, tmax_avF, tmax_stddevF, tempmax_avF, tempmax_stddevF = max_temp(LF, filenameF, quadraticnotcubicfit, toplotornot, usegaussian, 1, cutlower16, cutupper16)
print "Done with F"

maxhcs                 = [bm_avA, bm_avB, bm_avC, bm_avD, bm_avE]
maxhcs_stdv            = [bm_stddevA, bm_stddevB, bm_stddevC, bm_stddevD, bm_stddevE]
betaformaxes           = [tmax_avA, tmax_avB, tmax_avC, tmax_avD, tmax_avE]
betaformaxes_stdv      = [tmax_stddevA, tmax_stddevB, tmax_stddevC, tmax_stddevD, tmax_stddevE]
tempformaxes           = [tempmax_avA, tempmax_avB, tempmax_avC, tempmax_avD, tempmax_avE]
tempformaxes_stdv      = [tempmax_stddevA, tempmax_stddevB, tempmax_stddevC, tempmax_stddevD, tempmax_stddevE]
vec1dL                 = [1./LA, 1./LB, 1./LC, 1./LD, 1./LE]
vecL                   = [LA, LB, LC, LD, LE]

if bool16!=0:
    maxhcs.append(bm_avF)
    maxhcs_stdv.append(bm_stddevF)
    betaformaxes.append(tmax_avF)
    betaformaxes_stdv.append(tmax_stddevF)
    tempformaxes.append(tempmax_avF)
    tempformaxes_stdv.append(tempmax_stddevF)
    vec1dL.append(1./LF)
    vecL.append(LF)

maxhcs                 = array(maxhcs)
maxhcs_stdv            = array(maxhcs_stdv)
betaformaxes           = array(betaformaxes)
betaformaxes_stdv      = array(betaformaxes_stdv)
tempformaxes           = array(tempformaxes)
tempformaxes_stdv      = array(tempformaxes_stdv)
vec1dL                 = array(vec1dL)

# Print to file 
outfile = open("maxheatcap_Jensen.txt","w")
for i in range(len(maxhcs)):
    outfile.write('%.16f %.16f \n'%(maxhcs[i], maxhcs_stdv[i]))
outfile.close()

outfileT = open("temp_at_maxheatcap_Jensen_005.txt","w")
for i in range(len(tempformaxes)):
    outfileT.write('%.16f %.16f \n'%(tempformaxes[i], tempformaxes_stdv[i]))
outfileT.close()


max1dL_inplot = (max(vec1dL)*1.02)

# Regular plot
figure(figsize=(6,5))
errorbar(vecL, maxhcs, yerr=maxhcs_stdv, fmt=None, capsize=2)
hold('on')
title(r'Maxes of the specific heat', fontsize=14)
xlabel(r'$L$', fontsize=20)
ylabel(r'$c_{V}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, max(vecL), (min(maxhcs)-max(maxhcs_stdv))*0.99, (max(maxhcs)+max(maxhcs_stdv))*1.01])
show()

# Just because I am curious...
figure(figsize=(6,5))
errorbar(vec1dL, tempformaxes, yerr=tempformaxes_stdv, fmt=None, capsize=2)
hold('on')
title(r'Maxes of the specific heat', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$T_{max}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, max1dL_inplot, (min(tempformaxes)-max(tempformaxes_stdv))*0.999, (max(tempformaxes)+max(tempformaxes_stdv))*1.001])
show()

# log-log plot
figure(figsize=(6,5))
loglog(vecL, maxhcs)
hold('on')
loglog(vecL, maxhcs, 'o')
title(r'Maxes of the specific heat', fontsize=14)
xlabel(r'$L$', fontsize=20)
ylabel(r'$c_{V}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, max(vecL), (min(maxhcs)-max(maxhcs_stdv))*0.99, (max(maxhcs)+max(maxhcs_stdv))*1.01])
show()


# What I did for magnsusc
nL          = len(vecL)
logL        = zeros(nL)
logcv       = zeros(nL)
logcv_stdv  = zeros(nL)
logt        = zeros(nL)
logt_stdv   = zeros(nL)
cvpstdv     = zeros(nL)
cvmstdv     = zeros(nL)
tpstdv      = zeros(nL)
tmstdv      = zeros(nL)
for i in range(nL):
    logL[i]        = log(vecL[i])
    logcv[i]       = log(maxhcs[i])
    logcv_stdv[i]  = log(maxhcs_stdv[i])
    logt[i]        = log(tempformaxes[i])
    logt_stdv[i]   = log(tempformaxes_stdv[i])
    cvpstdv[i]     = maxhcs[i]+maxhcs_stdv[i]
    cvmstdv[i]     = maxhcs[i]-maxhcs_stdv[i]
    tpstdv[i]      = tempformaxes[i]+tempformaxes_stdv[i]
    tmstdv[i]      = tempformaxes[i]-tempformaxes_stdv[i]

comparewithlinex = array([vecL[0],vecL[nL-1]])
comparewithliney = array([maxhcs[0], maxhcs[nL-1]])

comparewithliney_pstdv = array([maxhcs[0]+maxhcs_stdv[0], maxhcs[nL-1]+maxhcs_stdv[nL-1]])
comparewithliney_mstdv = array([maxhcs[0]-maxhcs_stdv[0], maxhcs[nL-1]-maxhcs_stdv[nL-1]])
comparewithliney_combstdv1 = array([maxhcs[0]+maxhcs_stdv[0], maxhcs[nL-1]-maxhcs_stdv[nL-1]])
comparewithliney_combstdv2 = array([maxhcs[0]-maxhcs_stdv[0], maxhcs[nL-1]+maxhcs_stdv[nL-1]])

comparewithliney_t = array([maxhcs[0], maxhcs[nL-1]])



figure(figsize=(6,5))
loglog(vecL, maxhcs, basex=2, basey=2, label='Line between each point') # But... the error bars...
hold('on')
loglog(comparewithlinex, comparewithliney_t, 'r', basex=2, basey=2, label='Line between endpoints')
loglog(vecL, maxhcs, 'o', basex=2, basey=2, label='Points') #, basex=2 # <-- Can use commands such as this one to plot the line
title(r'$T$ at $\chi_{max}$', fontsize=14)
xlabel(r'$\log(L)$', fontsize=20)
ylabel(r'$\log(c_{V,max})$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower left")
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()

figure(figsize=(6,5))
loglog(vecL, maxhcs, basex=2, basey=2, label=r'T at $\chi_{max}$') # But... the error bars...
hold('on')
loglog(vecL, maxhcs, 'o', basex=2, basey=2)
loglog(vecL, cvpstdv, basex=2, basey=2, label=r'T+$\Delta\chi_{max}$ at $\chi_{max}$')
loglog(vecL, cvpstdv, 'o', basex=2, basey=2)
loglog(vecL, cvmstdv, basex=2, basey=2, label=r'T-$\Delta\chi_{max}$ at $\chi_{max}$')
loglog(vecL, cvmstdv, 'o', basex=2, basey=2) #, basex=2 # <-- Can use commands such as this one to plot the line
title(r'$T$ at $\chi_{max}$', fontsize=14)
xlabel(r'$\log(L)$', fontsize=20)
ylabel(r'$\log(c_{V,max})$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower left")
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()


# This line seems very straight, so we go ahead and use the endpoints of the log arrays:
print "From T at chi max:"
a,b = findline(logL[0], logL[nL-1], logcv[0], logcv[nL-1])
print "The exponent, gamma/nu: ", a
print "The weight: ", b

print " "
print "From T at chi max+stdv:"
a,b = findline(logL[0], logL[nL-1], log(maxhcs[0]+maxhcs_stdv[0]), log(maxhcs[nL-1]+maxhcs_stdv[nL-1]))
print "The exponent, gamma/nu: ", a
print "The weight: ", b

print " "
print "From T at chi max-stdv:"
a,b = findline(logL[0], logL[nL-1], log(maxhcs[0]-maxhcs_stdv[0]), log(maxhcs[nL-1]-maxhcs_stdv[nL-1]))
print "The exponent, gamma/nu: ", a
print "The weight: ", b

print " "
print "From T at chi max-stdv, max+stdv:"
a,b = findline(logL[0], logL[nL-1], log(maxhcs[0]-maxhcs_stdv[0]), log(maxhcs[nL-1]+maxhcs_stdv[nL-1]))
print "The exponent, gamma/nu: ", a
print "The weight: ", b

print " "
print "From T at chi max+stdv, max-stdv:"
a,b = findline(logL[0], logL[nL-1], log(maxhcs[0]+maxhcs_stdv[0]), log(maxhcs[nL-1]-maxhcs_stdv[nL-1]))
print "The exponent, gamma/nu: ", a
print "The weight: ", b
