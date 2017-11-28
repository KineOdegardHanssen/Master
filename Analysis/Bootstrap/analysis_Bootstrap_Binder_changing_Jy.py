from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys
from random import randint
from scipy.optimize import curve_fit

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
    
def giveexponential(a, b, c, x):
    lenx = len(x)
    expfunc = zeros(lenx)
    for i in range(lenx):
        expfunc[i] = a*exp(x[i]*b)+c
    return expfunc

def findintersection(beta, f1, f2):
    diff = f1 - f2
    for i in range(len(diff)):
        if (i+1)<len(diff): # Just a safeguard. We're doomed if this happens anyway
            if (diff[i] == 0) or (diff[i]*diff[i+1] < 0):
                return  beta[i], f1[i]
    #print "NB! Intersection not found!"
    return beta[0], f1[0]

def findline_importantlast(xvec, yvec):
    # Makes a line of the two last points in an array, to fit with
    length = len(yvec)
    a = (yvec[length-1]-yvec[length-2])/(xvec[length-1]-xvec[length-2]) # slope, a = (f(x2)-f(x1))/(x2-x1)
    b = yvec[length-1] - (a*xvec[length-1])                             # intersection, b = f(x)-ax
    return a, b

def givediff(xvec, yvec):
    # Creating a an array to compare to line fit 
    length = len(yvec)
    a = (yvec[length-1]-yvec[length-2])/(xvec[length-1]-xvec[length-2]) # slope, a = (f(x2)-f(x1))/(x2-x1)
    b = yvec[length-1] - (a*xvec[length-1])                             # intersection, b = f(x)-ax
    linepoint = a*xvec[length-3]+b
    diff = yvec[length-3]-linepoint
    reldiff = diff/yvec[length-3]
    diff2 = 0; reldiff2 = 0
    if length>3:
        linepoint2 = a*xvec[length-4]+b
        diff2 = yvec[length-4]-linepoint2
        reldiff2 = diff/yvec[length-4]
    rms = 0
    for i in range(length-1):
        rms += (yvec[i]-(a*xvec[i]+b))**2
    rms = sqrt(rms/length)
    return diff, reldiff, diff2, reldiff2, rms

def giverms(datavec, fitvec):
    length = len(datavec)
    rms = 0
    for i in range(length):
        rms += (fitvec[i]-datavec[i])**2
    rms = sqrt(rms/length)
    return rms

def findminima_abs(vec):
    lv = len(vec)
    minelem = 10000
    index = 10000
    for i in range(lv):
        #print "Element in vector:", vec[i]
        if (abs(vec[i])<minelem):
            minelem = abs(vec[i])
            index = i
            #print "In loop. Min elem:", minelem, "i:", i
    #print "Min elem:", minelem, "i:", index
    return minelem, index

def converttemp(crossingbeta):
    return 11.6045221/crossingbeta


def extract_crossing(LA, LB, filenameA, filenameB, quadraticnotcubicfit, toplotornot, cutit, cutlower, cutupper, deltaL): # Extract the beta value of the crossing between two graphs
    infileA = open(filenameA, "r")
    # Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
    
    # Sorting arrays
    betasA                 = []          # List of beta values
    binofbeta_firstindexA  = []          # I guess nested lists are a bit of a hassle
    binofbeta_lastindexA   = []          # Could also just feed the number of bins in. But, the less such input, the better, probably
    # Is the above unsafe? Should I just find the no of bins and add instead, a la, i*nbins+something?
    
    # Other quantities
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
    
    # Then B
    
    infileB = open(filenameB, "r")
    # Getting lists ready to store the data # Hmmm, now I have 3D spins...# Take the vector?
    
    # Sorting arrays
    betasB                 = []          # List of beta values
    binofbeta_firstindexB  = []          # I guess nested lists are a bit of a hassle
    binofbeta_lastindexB   = []          # Could also just feed the number of bins in. But, the less such input, the better, probably
    # Is the above unsafe? Should I just find the no of bins and add instead, a la, i*nbins+something?
    
    # Other quantities
    mxsc_bavsB        = []          # Index 4
    mxsc_abs_bavsB    = []          # Index 5
    mxssqc_bavsB      = []          # Index 6
    mxsquadc_bavsB    = []          # Index 7
    mysc_bavsB        = []          # Index 8
    mysc_abs_bavsB    = []          # Index 9
    myssqc_bavsB      = []          # Index 10
    mysquadc_bavsB    = []          # Index 11
    mzsc_bavsB        = []          # Index 12
    mzsc_abs_bavsB    = []          # Index 13
    mzssqc_bavsB      = []          # Index 14
    mzsquadc_bavsB    = []          # Index 15
    
    # Don't think I need the magnetizations of q=(0,0,0), but we'll see...
    # Read the rest of the lines
    lines = infileB.readlines()  # This works, I have checked it
    
    i = 0 # Some sort of counter to get the indices right
    betabefore = 0 # Test the betas so we can gather the bins.
    binofbeta_firstindexB.append(0) 
    # Should I have some list over the indices of the betas belonging to each bin?
    if toplotornot==0:
        plotter = 1
    else:
        plotter = 0
    
    # Getting data from the file
    for line in lines:
        words = line.split()
        if len(words) != 0:
            # Betas
            beta = float(words[0])
            if beta != betabefore:
                betasB.append(beta)
                betabefore = beta
                if i != 0:
                    binofbeta_firstindexB.append(i)
                    binofbeta_lastindexB.append(i-1)
            # <m_x(q)>
            en = float(words[4])
            mxsc_bavsB.append(en)            # Index 4
            # <|m_x(q)|>
            en = float(words[5])
            mxsc_abs_bavsB.append(en)        # Index 5
            # <m^2_x(q)>
            en = float(words[6])
            mxssqc_bavsB.append(en)          # Index 6
            # <m^4_x(q)>
            en = float(words[7])
            mxsquadc_bavsB.append(en)        # Index 7
            # <m_y(q)>
            en = float(words[8])
            mysc_bavsB.append(en)            # Index 8
            # <|m_y(q)|>
            en = float(words[9])
            mysc_abs_bavsB.append(en)        # Index 9
            # <m^2_y(q)>
            en = float(words[10])
            myssqc_bavsB.append(en)          # Index 10
            # <m^4_y(q)>
            en = float(words[11])
            mysquadc_bavsB.append(en)        # Index 11
            # <m_z(q)>
            en = float(words[12])
            mzsc_bavsB.append(en)            # Index 12
            # <|m_z(q)|>
            en = float(words[13])
            mzsc_abs_bavsB.append(en)        # Index 13
            # <m^2_z(q)>
            en = float(words[14])
            mzssqc_bavsB.append(en)          # Index 14
            # <m^4_z(q)>
            en = float(words[15])
            mzsquadc_bavsB.append(en)        # Index 15
        i += 1 # Increasing the counter
        
    binofbeta_lastindexB.append(i-1)  # The index has been updated one time too many
    
    # Remember to close the file
    infileB.close()
    
    nbinsB = zeros(len(betasB))
    for i in range(0,len(betasB)):
        nbinsB[i] = binofbeta_lastindexB[i]-binofbeta_firstindexB[i]+1
    
    no_of_bins_each_betaB = binofbeta_lastindexB[0]+1 # Due to the indexing in Python + we don't allow
    no_of_betasB          = len(betasB)               # The number of different temperatures ### Should I convert to T now?
    no_of_BootstraprunsB  = no_of_BootstraprunsA      # Probably OK this way... Good to be flexible...
    
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
            if (betasA[i]<=betamax and betasA[i]>=betamin):
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
        # For B # Well, A and B must be of the same size, I think, but this can't hurt
        smallbetaarray = []
        binofbeta_firstindex_thosewewant = []
        binofbeta_lastindex_thosewewant = []
        no_of_betas_wewant= 0
        for i in range(0, no_of_betasB):
            if (betasB[i]<=betamax and betasB[i]>=betamin):
                smallbetaarray.append(betasB[i])
                nbins_wewant.append(nbinsB[i])
                binofbeta_firstindex_thosewewant.append(binofbeta_firstindexB[i])
                binofbeta_lastindex_thosewewant.append(binofbeta_lastindexB[i])
                no_of_betas_wewant+=1
        no_of_betasB = no_of_betas_wewant
        betasB = smallbetaarray
        nbinsB = nbins_wewant
        binofbeta_firstindexB = binofbeta_firstindex_thosewewant
        binofbeta_lastindexB = binofbeta_lastindex_thosewewant
    betasA = array(betasA) # Not sure if this is neccessary, but maybe...
    
    intersectionbeta_av = 0
    intersectionbetas = zeros(no_of_BootstraprunsA)
    intersectiontemp_av = 0
    intersectiontemps = zeros(no_of_BootstraprunsA)
    intersectionU_av = 0
    intersectionUs = zeros(no_of_BootstraprunsA)

    #print "L=", LA, " and ", LB
    #print "no. of betas, A: ", no_of_betasA
    #print "no. of betas, B: ", no_of_betasB
    for k in range(0, no_of_BootstraprunsA): # Want to repeat the random selection of bins a number of times
        #bindervalueszcompA = bootstrapping_each(no_of_betasA, binofbeta_firstindexA, binofbeta_lastindexA, nbinsA, mxsc_bavsA, mysc_bavsA, mzsc_bavsA, mxsc_abs_bavsA, mysc_abs_bavsA, mzsc_abs_bavsA, mxssqc_bavsA, myssqc_bavsA, mzssqc_bavsA, mxsquadc_bavsA, mysquadc_bavsA, mzsquadc_bavsA)
        #bindervalueszcompB = bootstrapping_each(no_of_betasB, binofbeta_firstindexB, binofbeta_lastindexB, nbinsB, mxsc_bavsB, mysc_bavsB, mzsc_bavsB, mxsc_abs_bavsB, mysc_abs_bavsB, mzsc_abs_bavsB, mxssqc_bavsB, myssqc_bavsB, mzssqc_bavsB, mxsquadc_bavsB, mysquadc_bavsB, mzsquadc_bavsB)
        
        bindervalueszcompA = []
        for i in range(0, no_of_betasA): # For doing it for every temperature
            ### Reset all quantities I want to find by bootstrap        
            # From A
            mxA = 0; mxabsA = 0; mx2A = 0; mx4A = 0;
            myA = 0; myabsA = 0; my2A = 0; my4A = 0;
            mzA = 0; mzabsA = 0; mz2A = 0; mz4A = 0;
            nbinsthis = int(nbinsA[i])
            for j in range(0, nbinsthis): # For finding the averages
                ###Draw a random integer in [binofbeta_firstindexA[i], binofbeta_lastindexA[i]].
                n = randint(binofbeta_firstindexA[i], binofbeta_lastindexA[i])  # Draw a random number n times
                ###Extract O(n), add to average function
                # From A
                mxA += mxsc_bavsA[n]; mxabsA += mxsc_abs_bavsA[n]; mx2A += mxssqc_bavsA[n]; mx4A += mxsquadc_bavsA[n];
                myA += mysc_bavsA[n]; myabsA += mysc_abs_bavsA[n]; my2A += myssqc_bavsA[n]; my4A += mysquadc_bavsA[n];
                mzA += mzsc_bavsA[n]; mzabsA += mzsc_abs_bavsA[n]; mz2A += mzssqc_bavsA[n]; mz4A += mzsquadc_bavsA[n];
            ### Divide the average of O(n) by no_of_bins_each_beta to ACTUALLY make the average
            # From A
            mxA = mxA/nbinsthis; mxabsA = mxabsA/nbinsthis; mx2A = mx2A/nbinsthis; mx4A = mx4A/nbinsthis;
            myA = myA/nbinsthis; myabsA = myabsA/nbinsthis; my2A = my2A/nbinsthis; my4A = my4A/nbinsthis;
            mzA = mzA/nbinsthis; mzabsA = mzabsA/nbinsthis; mz2A = mz2A/nbinsthis; mz4A = mz4A/nbinsthis;
            ### Find the Binder cumulant
            # For A
            thisbinderzA = 1 - (mz4A/(3.*mz2A*mz2A))
            bindervalueszcompA.append(thisbinderzA)
        
        bindervalueszcompB = []
        for i in range(0, no_of_betasB): # For doing it for every temperature
            ### Reset all quantities I want to find by bootstrap        
            # From B
            mxB = 0; mxabsB = 0; mx2B = 0; mx4B = 0;
            myB = 0; myabsB = 0; my2B = 0; my4B = 0;
            mzB = 0; mzabsB = 0; mz2B = 0; mz4B = 0;
            nbinsthis = int(nbinsB[i])
            for j in range(0, nbinsthis): # For finding the averages
                ###Draw a random integer in [binofbeta_firstindexA[i], binofbeta_lastindexA[i]].
                n = randint(binofbeta_firstindexB[i], binofbeta_lastindexB[i])  # Draw a random number n times
                ###Extract O(n), add to average function
                # From B
                mxB += mxsc_bavsB[n]; mxabsB += mxsc_abs_bavsB[n]; mx2B += mxssqc_bavsB[n]; mx4B += mxsquadc_bavsB[n];
                myB += mysc_bavsB[n]; myabsB += mysc_abs_bavsB[n]; my2B += myssqc_bavsB[n]; my4B += mysquadc_bavsB[n];
                mzB += mzsc_bavsB[n]; mzabsB += mzsc_abs_bavsB[n]; mz2B += mzssqc_bavsB[n]; mz4B += mzsquadc_bavsB[n];
            ### Divide the average of O(n) by no_of_bins_each_beta to ACTUALLY make the average
            # From B
            mxB = mxB/nbinsthis; mxabsB = mxabsB/nbinsthis; mx2B = mx2B/nbinsthis; mx4B = mx4B/nbinsthis;
            myB = myB/nbinsthis; myabsB = myabsB/nbinsthis; my2B = my2B/nbinsthis; my4B = my4B/nbinsthis;
            mzB = mzB/nbinsthis; mzabsB = mzabsB/nbinsthis; mz2B = mz2B/nbinsthis; mz4B = mz4B/nbinsthis;
            ### Find the Binder cumulant
            # For B
            thisbinderzB = 1 - (mz4B/(3.*mz2B*mz2B))
            bindervalueszcompB.append(thisbinderzB)
            #if plotter==0:
            #    print "L=", LA, " and ", LB, ": bindervalueszcompB[",i,"]=", bindervalueszcompB[i], " which should be ", thisbinderzB
        
        ### Do a fitting of the results. To a line, a quadr. func., a pol. of the third degree. Check this. OFS suggested quadr
        if quadraticnotcubicfit==0:
            # For A
            bindervalueszcompA = array(bindervalueszcompA) 
            fitvectorzcompA = polyfit(betasA, bindervalueszcompA, 1) # Fits the function points to a line
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1];
            # For B
            bindervalueszcompB = array(bindervalueszcompB) 
            fitvectorzcompB = polyfit(betasB, bindervalueszcompB, 1) # Fits the function points to a line
            azB = fitvectorzcompB[0]; bzB = fitvectorzcompB[1]; 
                        
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = giveline(azA, bzA, fbetas) # Should this really be quadratic? I look at a small interval...
            fB = giveline(azB, bzB, fbetas)
        if quadraticnotcubicfit==1:
            # For A
            bindervalueszcompA = array(bindervalueszcompA) 
            fitvectorzcompA = polyfit(betasA, bindervalueszcompA, 2) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]
            # For B
            bindervalueszcompB = array(bindervalueszcompB) 
            fitvectorzcompB = polyfit(betasB, bindervalueszcompB, 2) # Fits the function points to a quadratic polynomial
            azB = fitvectorzcompB[0]; bzB = fitvectorzcompB[1]; czB = fitvectorzcompB[2]
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givequadratic(azA, bzA, czA, fbetas)
            fB = givequadratic(azB, bzB, czB, fbetas)
        if quadraticnotcubicfit==2:
            # For A
            bindervalueszcompA = array(bindervalueszcompA) 
            fitvectorzcompA = polyfit(betasA, bindervalueszcompA, 3) # Fits the function points to a quadratic polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]; dzA = fitvectorzcompA[3];
            # For B
            bindervalueszcompB = array(bindervalueszcompB) 
            fitvectorzcompB = polyfit(betasB, bindervalueszcompB, 3) # Fits the function points to a quadratic polynomial
            azB = fitvectorzcompB[0]; bzB = fitvectorzcompB[1]; czB = fitvectorzcompB[2]; dzB = fitvectorzcompB[3];
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givecubic(azA, bzA, czA, dzA, fbetas)
            fB = givecubic(azB, bzB, czB, dzB, fbetas)
        if quadraticnotcubicfit==3:
            # For A
            bindervalueszcompA = array(bindervalueszcompA) 
            fitvectorzcompA = polyfit(betasA, bindervalueszcompA, 4) # Fits the function points to a quadruple polynomial
            azA = fitvectorzcompA[0]; bzA = fitvectorzcompA[1]; czA = fitvectorzcompA[2]; dzA = fitvectorzcompA[3];ezA = fitvectorzcompA[4];
            # For B
            bindervalueszcompB = array(bindervalueszcompB) 
            fitvectorzcompB = polyfit(betasB, bindervalueszcompB, 4) # Fits the function points to a quadratic polynomial
            azB = fitvectorzcompB[0]; bzB = fitvectorzcompB[1]; czB = fitvectorzcompB[2]; dzB = fitvectorzcompB[3]; ezB = fitvectorzcompB[4];
            
            nfbetas = 1000
            fbetas = linspace(betasA[0], betasA[len(betasA)-1], nfbetas)
            
            fA = givequadruple(azA, bzA, czA, dzA, ezA, fbetas) # Should this really be quadruple? I look at a small interval...
            fB = givequadruple(azB, bzB, czB, dzB, ezB, fbetas)
        crossingbeta, crossingU = findintersection(fbetas, fA, fB)
        intersectionbeta_av += crossingbeta
        intersectionbetas[k] = crossingbeta
        temp = converttemp(crossingbeta)
        intersectiontemp_av += temp
        intersectiontemps[k] = temp
        intersectionU_av += crossingU
        intersectionUs[k] = crossingU
        
        if plotter==0:
            figure(figsize=(6,5))
            plot(fbetas, fA, label='L=%i, fit'%LA)
            hold('on')
            plot(fbetas, fB, label='L=%i, fit'%LB)
            plot(betasA, bindervalueszcompA, 'ro', label='L=%i, data'%LA)
            plot(betasB, bindervalueszcompB, 'bo', label='L=%i, data'%LB)
            title(r'Crossings between Binder cumulant graphs L=%i and %i'%(LA, LB), fontsize=14)
            xlabel(r'$\beta$', fontsize=20)
            ylabel(r'$U_{L,z}$', fontsize=20)
            legend(loc="upper left")
            tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
            show()
            print "L=", LA, "and", LB, ", beta = ", crossingbeta
            if k==2: # Getting more plots
                plotter = 1
        
    intersectionbeta_av = intersectionbeta_av/no_of_BootstraprunsA
    intersectionbeta_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        intersectionbeta_stddev += (intersectionbeta_av-intersectionbetas[i])*(intersectionbeta_av-intersectionbetas[i])
    intersectionbeta_stddev = sqrt(intersectionbeta_stddev/(no_of_BootstraprunsA-1))
    
    intersectiontemp_av = intersectiontemp_av/no_of_BootstraprunsA
    intersectiontemp_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        intersectiontemp_stddev += (intersectiontemp_av-intersectiontemps[i])*(intersectiontemp_av-intersectiontemps[i])
    intersectiontemp_stddev = sqrt(intersectiontemp_stddev/(no_of_BootstraprunsA-1))
    
    intersectionU_av = intersectionU_av/no_of_BootstraprunsA
    intersectionU_stddev = 0
    for i in range(0, no_of_BootstraprunsA):
        intersectionU_stddev += (intersectionU_av-intersectionUs[i])*(intersectionU_av-intersectionUs[i])
    intersectionU_stddev = sqrt(intersectionU_stddev/(no_of_BootstraprunsA-1))
    
    return intersectionbeta_av, intersectionbeta_stddev, intersectiontemp_av, intersectiontemp_stddev, intersectionU_av, intersectionU_stddev

# Files, L = 12
L12fileJy0p000 = "fcc12x12x12yopen_closetophasetransition_nnnJy0p0_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs_big.txt"
L12fileJy0p300 = "fcc12x12x12yopen_closetophasetransition_nnnJy0p3_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs_big.txt"
L12fileJy0p500 = "fcc12x12x12yopen_closetophasetransition_nnnJy0p5_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs_big.txt"
L12fileJy0p600 = "fcc12x12x12yopen_closetophasetransition_nnnJy0p6_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs_big.txt" #
L12fileJy0p660 = "fcc12x12x12yopen_TvsJy_phasediagr_Jy0p66_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs.txt"
L12fileJy0p670 = "fcc12x12x12yopen_ds_shortrange_beta0p782to0p792_Nbeta10_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_slowcool_binavgs.txt"
L12fileJy0p671 = "fcc12x12x12yopen_TvsJy_phasediagr_Jy0p671_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_divseeds_slowcool_binavgs.txt"
L12fileJy0p675 = "fcc12x12x12yopen_TvsJy_phasediagr_Jy0p675_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_divseeds_slowcool_binavgs.txt"
L12fileJy0p700 = "fcc12x12x12yopen_TvsJy_phasediagr_Jy0p7_SMALLCAT_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_divseeds_slowcool_binavgs.txt"
#L12fileJy0p700 = "fcc12x12x12yopen_TvsJy_phasediagr_Jy0p70_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs.txt" # This is bad. Use the one above


# Files, L = 14
L14fileJy0p000 = "fcc14x14x14yopen_closetophasetransition_nnnJy0p0_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs_big.txt"
L14fileJy0p300 = "fcc14x14x14yopen_closetophasetransition_nnnJy0p3_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs_big.txt"
L14fileJy0p500 = "fcc14x14x14yopen_closetophasetransition_nnnJy0p5_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs_big.txt"
L14fileJy0p600 = "fcc14x14x14yopen_closetophasetransition_nnnJy0p6_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs_big.txt" #
L14fileJy0p660 = "fcc14x14x14yopen_TvsJy_phasediagr_Jy0p66_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs.txt"
L14fileJy0p670 = "fcc14x14x14yopen_ds_shortrange_beta0p782to0p792_Nbeta10_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_slowcool_binavgs.txt"
L14fileJy0p671 = "fcc14x14x14yopen_TvsJy_phasediagr_Jy0p671_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_divseeds_slowcool_binavgs.txt"
L14fileJy0p675 = "fcc14x14x14yopen_TvsJy_phasediagr_Jy0p675_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_divseeds_slowcool_binavgs.txt"
L14fileJy0p700 = "fcc14x14x14yopen_TvsJy_phasediagr_Jy0p7_SMALLCAT_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc10000_bins100_divseeds_slowcool_binavgs.txt"
#L14fileJy0p700 = "fcc14x14x14yopen_TvsJy_phasediagr_Jy0p70_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc100000_bins100_divseeds_slowcool_binavgs.txt" # This is bad. Use the one above
# Bools for degree of fit, plotting, range of values
quadraticnotcubicfit = 2 # 0 gives a line fit to the Binder cumulants, 1 gives a quadratic fit, 2 gives cubic, 3 gives quadruple
toplotornot = 0 # 0, don't plot, anything else leads to a plot

# If we only want to look at a small temperature interval
cutit = 0 # 0: Find a fit for the whole beta range; 1: Find a fit for a small interval

cutlowerAB = 0.777 # Just to pass an argument
cutupperAB = 0.784

cutlower0p0 = 0.553
cutupper0p0 = 0.567

cutlower0p3 = 0.62
cutupper0p3 = 0.64

cutlower0p5 = 0.692
cutupper0p5 = 0.708

cutlower0p6 = 0.737
cutupper0p6 = 0.754

cutlower0p66 = 0.777
cutupper0p66 = 0.784

cutlower0p67 = 0.782 # Already cut
cutupper0p67 = 0.792

cutlower0p671 = 0.784
cutupper0p671 = 0.792

cutlower0p675 = 0.784
cutupper0p675 = 0.796

#cutlower1 = cutlower2
#cutupper1 = cutupper2

deltaL = 2 # Difference in the size of the systems for the graphs we use for the intersection
LA = 12
LB = 14

print "starting Bootstrap, Jy=0:"
ib_av0p000, ib_stdv0p000, iT_av0p000, iT_stdv0p000, iU_av0p000, iU_stdv0p000 = extract_crossing(LA, LB, L12fileJy0p000, L14fileJy0p000, quadraticnotcubicfit, toplotornot, 1, cutlower0p0, cutupper0p0, deltaL)
print "starting Bootstrap, Jy=0.3:"
ib_av0p300, ib_stdv0p300, iT_av0p300, iT_stdv0p300, iU_av0p300, iU_stdv0p300 = extract_crossing(LA, LB, L12fileJy0p300, L14fileJy0p300, quadraticnotcubicfit, toplotornot, 1, cutlower0p3, cutupper0p3, deltaL)
print "starting Bootstrap, Jy=0.5:"
ib_av0p500, ib_stdv0p500, iT_av0p500, iT_stdv0p500, iU_av0p500, iU_stdv0p500 = extract_crossing(LA, LB, L12fileJy0p500, L14fileJy0p500, quadraticnotcubicfit, toplotornot, 1, cutlower0p5, cutupper0p5, deltaL)
print "starting Bootstrap, Jy=0.6:"
ib_av0p600, ib_stdv0p600, iT_av0p600, iT_stdv0p600, iU_av0p600, iU_stdv0p600 = extract_crossing(LA, LB, L12fileJy0p600, L14fileJy0p600, quadraticnotcubicfit, toplotornot, 1, cutlower0p6, cutupper0p6, deltaL)
print "starting Bootstrap, Jy=0.66:"
ib_av0p660, ib_stdv0p660, iT_av0p660, iT_stdv0p660, iU_av0p660, iU_stdv0p660 = extract_crossing(LA, LB, L12fileJy0p660, L14fileJy0p660, quadraticnotcubicfit, toplotornot, 1, cutlower0p66, cutupper0p66, deltaL)
print "starting Bootstrap, Jy=0.67:"
ib_av0p670, ib_stdv0p670, iT_av0p670, iT_stdv0p670, iU_av0p670, iU_stdv0p670 = extract_crossing(LA, LB, L12fileJy0p670, L14fileJy0p670, quadraticnotcubicfit, toplotornot, cutit, cutlowerAB, cutupperAB, deltaL)
print "starting Bootstrap, Jy=0.:671"
ib_av0p671, ib_stdv0p671, iT_av0p671, iT_stdv0p671, iU_av0p671, iU_stdv0p671 = extract_crossing(LA, LB, L12fileJy0p671, L14fileJy0p671, quadraticnotcubicfit, toplotornot, 1, cutlower0p671, cutupper0p671, deltaL)
print "starting Bootstrap, Jy=0.675:"
ib_av0p675, ib_stdv0p675, iT_av0p675, iT_stdv0p675, iU_av0p675, iU_stdv0p675 = extract_crossing(LA, LB, L12fileJy0p675, L14fileJy0p675, quadraticnotcubicfit, toplotornot, 1, cutlower0p675, cutupper0p675, deltaL)
#ib_av0p700, ib_stdv0p700, iT_av0p700, iT_stdv0p700, iU_av0p700, iU_stdv0p700 = extract_crossing(LA, LB, L12fileJy0p700, L14fileJy0p700, quadraticnotcubicfit, toplotornot, cutit, cutlowerAB, cutupperAB, deltaL)


Jys = array([0,0.3,0.5, 0.6, 0.66, 0.67, 0.671, 0.675])#, 0.7

temps_av = array([iT_av0p000, iT_av0p300, iT_av0p500, iT_av0p600, iT_av0p660, iT_av0p670, iT_av0p671, iT_av0p675]) #, iT_av0p700

temps_stdv = array([iT_stdv0p000,iT_stdv0p300,  iT_stdv0p500, iT_stdv0p600, iT_stdv0p660, iT_stdv0p670, iT_stdv0p671,  iT_stdv0p675]) #, iT_stdv0p700

 
Us_av = array([iU_av0p000,iU_av0p300,iU_av0p500,iU_av0p600,iU_av0p660,iU_av0p670,iU_av0p671,iU_av0p675])
Us_stdv = array([iU_stdv0p000,iU_stdv0p300,iU_stdv0p500,iU_stdv0p600,iU_stdv0p660,iU_stdv0p670,iU_stdv0p671,iU_stdv0p675])


# A closer look
closeup_Jys = array([0.66, 0.67, 0.671, 0.675])
closeup_t      = array([iT_av0p660, iT_av0p670, iT_av0p671, iT_av0p675])
closeup_t_stdv = array([iT_stdv0p660, iT_stdv0p670, iT_stdv0p671, iT_stdv0p675])
closeup_U      = array([iU_av0p660,iU_av0p670,iU_av0p671,iU_av0p675])
closeup_U_stdv = array([iU_stdv0p660,iU_stdv0p670,iU_stdv0p671,iU_stdv0p675])

lent = len(Jys)
outfilecrt = open("BinderPT_T_001.txt", 'w')
outfilecrU = open("BinderPT_U_001.txt", 'w')
for i in range(lent):
    outfilecrt.write('%.16f %.16f \n'%(temps_av[i], temps_stdv[i])) # Guess I don't need the standard deviations...
    outfilecrU.write('%.16f %.16f \n'%(Us_av[i], Us_stdv[i])) # Guess I don't need the standard deviations...
outfilecrt.close()
outfilecrU.close()

lent = len(closeup_Jys)
outfilecrt2 = open("BinderPT_T_closeup_001.txt", 'w')
outfilecrU2 = open("BinderPT_U_closeup_001.txt", 'w')
for i in range(lent):
    outfilecrt2.write('%.16f %.16f \n'%(closeup_t[i], closeup_t_stdv[i])) # Guess I don't need the standard deviations...
    outfilecrU2.write('%.16f %.16f \n'%(closeup_U[i], closeup_U_stdv[i])) # Guess I don't need the standard deviations...
outfilecrt2.close()
outfilecrU2.close()

figure(figsize=(6,5))
errorbar(Jys, temps_av, yerr=temps_stdv, fmt=None, capsize=2)
hold('on')
title(r'Binder crossings vs $J_y$, L=12 and L=14 ', fontsize=14)
xlabel(r'$J_y$', fontsize=20)
ylabel(r'T', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([-0.05, 0.7, 14, 21])
show()

figure(figsize=(6,5))
errorbar(Jys, temps_av, yerr=temps_stdv, capsize=2)
hold('on')
#plot(Jys, temps_av, 'o')
title(r'Binder crossings vs $J_y$, L=12 and L=14 ', fontsize=14)
xlabel(r'$J_y$', fontsize=20)
ylabel(r'T', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([-0.05, 0.7, 14, 21])
show()

figure(figsize=(6,5))
errorbar(Jys, Us_av, yerr=Us_stdv, fmt=None, capsize=2)
hold('on')
title(r'Binder crossings vs $J_y$, L=12 and L=14 ', fontsize=14)
xlabel(r'$J_y$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([-0.05, 0.73, (min(Us_av)-max(Us_stdv))*0.999, (max(Us_av)+max(Us_stdv))*1.001])
show()

figure(figsize=(6,5))
errorbar(Jys, Us_av, yerr=Us_stdv, capsize=2)
hold('on')
#plot(Jys, Us_av, 'o')
title(r'Binder crossings vs $J_y$, L=12 and L=14 ', fontsize=14)
xlabel(r'$J_y$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([-0.05, 0.73, (min(Us_av)-max(Us_stdv))*0.999, (max(Us_av)+max(Us_stdv))*1.001])
show()

###### Closeup

figure(figsize=(6,5))
errorbar(closeup_Jys, closeup_t, yerr=closeup_t_stdv, fmt=None, capsize=2)
hold('on')
title(r'Binder crossings vs $J_y$, L=12 and L=14 ', fontsize=14)
xlabel(r'$J_y$', fontsize=20)
ylabel(r'T', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0.659, 0.676, 14.65, 14.9])
show()

figure(figsize=(6,5))
errorbar(closeup_Jys, closeup_t, yerr=closeup_t_stdv, capsize=2)
hold('on')
#plot(closeup_Jys, closeup_t, 'o')
title(r'Binder crossings vs $J_y$, L=12 and L=14 ', fontsize=14)
xlabel(r'$J_y$', fontsize=20)
ylabel(r'T', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0.659, 0.676, 14.65, 14.9])
show()

figure(figsize=(6,5))
errorbar(closeup_Jys, closeup_U, yerr=closeup_U_stdv, fmt=None, capsize=2)
hold('on')
title(r'Binder crossings vs $J_y$, L=12 and L=14 ', fontsize=14)
xlabel(r'$J_y$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0.659, 0.676, 0.4, 0.414])
show()

figure(figsize=(6,5))
errorbar(closeup_Jys, closeup_U, yerr=closeup_U_stdv, capsize=2)
hold('on')
#plot(closeup_Jys, closeup_U, 'o')
title(r'Binder crossings vs $J_y$, L=12 and L=14 ', fontsize=14)
xlabel(r'$J_y$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0.659, 0.676, 0.4, 0.414])
show()
