from numpy import *
from scipy import *
from matplotlib.pyplot import *
import time
import sys
from random import randint
from scipy.optimize import curve_fit

# I will extract only info about Tc from the Binder cumulant first

def giveline(a,b,x):
    outpol = zeros(len(x)) 
    for i in range(len(x)):
        outpol[i] = a*x[i]+b
    return outpol

def givequadratic(a, b, c, x):
    outpol = zeros(len(x)) 
    for i in range(len(x)):
        outpol[i] = a*x[i]**2+b*x[i]+c
    return outpol
    
def givecubic(a,b,c,d,x):
    outpol = zeros(len(x)) 
    for i in range(len(x)):
        outpol[i] = a*x[i]**3+b*x[i]**2+c*x[i]+d
    return outpol

def givequadruple(a,b,c,d,e,x):
    outpol = zeros(len(x)) 
    for i in range(len(x)):
        outpol[i] = a*x[i]**4+b*x[i]**3+c*x[i]**2+d*x[i]+e
    return outpol
    
def givepowerlaw(k,n,x):
    return k*x**n
    
def giveexponential(a, b, c, x):
    lenx = len(x)
    expfunc = zeros(lenx)
    for i in range(lenx):
        expfunc[i] = a*exp(x[i]*b)+c
    return expfunc

def fittingfunc(x, crossing, weight, exponent): # x = L**(-1)
    return crossing + weight*x**(exponent)

def give_fittingfunc(x, crossing, weight, exponent): # x = L**(-1)
    lenx = len(x)
    expfunc = zeros(lenx)
    for i in range(lenx):
        expfunc[i] = crossing + weight*x[i]**(exponent)
    return expfunc

def givevalueofderivative(point, coefficients):
    no_c = len(coefficients) # Number of coefficients
    value_derivative = 0
    for i in range(no_c-1):
        value_derivative = (i+1)*(point**i)*coefficients[i+1]
    
def find_1dnu_L(L, deltaL, crossingbeta, crossingbeta_stdv, coeffsA, coeffs_stdA, coeffsB, coeffs_stdB):
    # L is LA
    # Setting a lot of quantities to make the equations more readable
    x = crossingbeta; dx = crossingbeta_stdv
    a1A = coeffsA[1]; a2A = coeffsA[2]; a3A = coeffsA[3]; a4A = coeffsA[4]
    a1B = coeffsB[1]; a2B = coeffsB[2]; a3B = coeffsB[3]; a4B = coeffsB[4]
    da1A = coeffs_stdA[1]; da2A = coeffs_stdA[2]; da3A = coeffs_stdA[3]; da4A = coeffs_stdA[4]
    da1B = coeffs_stdB[1]; da2B = coeffs_stdB[2]; da3B = coeffs_stdB[3]; da4B = coeffs_stdB[4]
    # Abbreviating polynomials
    fA = a1A + 2*a2A*x + 3*a3A*x**2 + 4*a4A*x**3
    fB = a1B + 2*a2B*x + 3*a3B*x**2 + 4*a4B*x**3
    dfA = 2*a2A + 6*a3A*x + 12*a4A*x**2
    dfB = 2*a2B + 6*a3B*x + 12*a4B*x**2
    print "fA = ", fA
    print "fB = ", fB
    print "dfA = ", dfA
    print "dfB = ", dfB
    # Setting another quanity I am going to use in both equations
    deltaLdL = float(deltaL)/L
    log1pdeltaLdL = log(1+deltaLdL)
    print "log(1+dL/L)= ",log1pdeltaLdL
    # Finding the estimate of the inverse exponent
    inverseexponent = log(abs(fB/fA))/log1pdeltaLdL       ####### This will be returned
    print "1/nu = ", inverseexponent
    # Finding the error estimate:
    # 1. Contribution from the uncertainty in the coefficients
    variance = (da1A/fA)**2+(2*x*da2A/fA)**2+(3*x**2*da3A/fA)**2+(4*x**3*da4A/fA)**2 + (da1B/fB)**2+(2*x*da2B/fB)**2+(3*x**2*da3B/fB)**2+(4*x**3*da4B/fB)**2
    # 2. Contribution from the uncertainty in the crossing temperature
    variance += ((dfB/fB - dfA/fA)*dx)**2
    # 3. Taking the square root and multiplying by the common factor of all the errors
    iexpstddev = sqrt(variance)/abs(log1pdeltaLdL)    ####### This will be returned
    print "delta(1/nu) = ", inverseexponent
    return inverseexponent, iexpstddev

def findintersection(beta, f1, f2):
    diff = f1 - f2
    for i in range(len(diff)):
        if (i+1)<len(diff): # Just a safeguard. We're doomed if this happens anyway
            if (diff[i] == 0) or (diff[i]*diff[i+1] < 0):
                return  beta[i], f1[i]
    print "NB! Intersection not found!"
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

def bootstrapping_each(no_of_betasA, binofbeta_firstindexA, binofbeta_lastindexA, nbinsA, mxsc_bavsA, mysc_bavsA, mzsc_bavsA, mxsc_abs_bavsA, mysc_abs_bavsA, mzsc_abs_bavsA, mxssqc_bavsA, myssqc_bavsA, mzssqc_bavsA, mxsquadc_bavsA, mysquadc_bavsA, mzsquadc_bavsA):
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
    return bindervalueszcompA


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
                #print "betaA = ", beta
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
        #print "Line no:", i
        if len(words) != 0:
            # Betas
            beta = float(words[0])
            if beta != betabefore:
                #print "betaB = ", beta
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

    #print "no_of_betasA =", no_of_betasA
    #print "no_of_betasB =", no_of_betasB
    
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
                #print "In cutting procedure. betasA[i] = ", betasA[i]
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
            #print "betasB[",i,"] =", betasB[i]
            if (betasB[i]<=betamax and betasB[i]>=betamin):
                #print "In cutting procedure. betasB[i] = ", betasB[i]
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

    av_a0A = 0
    av_a1A = 0
    av_a2A = 0
    av_a3A = 0
    av_a4A = 0
    each_a0A = zeros(no_of_BootstraprunsA)
    each_a1A = zeros(no_of_BootstraprunsA)
    each_a2A = zeros(no_of_BootstraprunsA)
    each_a3A = zeros(no_of_BootstraprunsA)
    each_a4A = zeros(no_of_BootstraprunsA)
    
    av_a0B = 0
    av_a1B = 0
    av_a2B = 0
    av_a3B = 0
    av_a4B = 0
    each_a0B = zeros(no_of_BootstraprunsA)
    each_a1B = zeros(no_of_BootstraprunsA)
    each_a2B = zeros(no_of_BootstraprunsA)
    each_a3B = zeros(no_of_BootstraprunsA)
    each_a4B = zeros(no_of_BootstraprunsA)

    #print "betasA = ", betasA

    exponent1dnu_av = 0
    exponent1dnus = zeros(no_of_BootstraprunsA)
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

            # Accumulating the average and storing the value of the coefficients
            av_a0A += azA
            av_a1A += bzA
            av_a0B += azB
            av_a1B += bzB
            each_a0A[k] = azA
            each_a1A[k] = bzA
            each_a0B[k] = azB
            each_a1B[k] = bzB
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

            # Accumulating the average and storing the value of the coefficients
            av_a0A += azA
            av_a1A += bzA
            av_a2A += czA
            av_a0B += azB
            av_a1B += bzB
            av_a2B += czB
            each_a0A[k] = azA
            each_a1A[k] = bzA
            each_a2A[k] = czA
            each_a0B[k] = azB
            each_a1B[k] = bzB
            each_a2B[k] = czB
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
            #print "len(fbetas)=", len(fbetas)
            #print fbetas
            fA = givecubic(azA, bzA, czA, dzA, fbetas)
            fB = givecubic(azB, bzB, czB, dzB, fbetas)

            #if LA==6:
            #    print "a0A = ", azA

            # Accumulating the average and storing the value of the coefficients
            av_a0A += azA
            av_a1A += bzA
            av_a2A += czA
            av_a3A += dzA
            av_a0B += azB
            av_a1B += bzB
            av_a2B += czB
            av_a3B += dzB
            each_a0A[k] = azA
            each_a1A[k] = bzA
            each_a2A[k] = czA
            each_a3A[k] = dzA
            each_a0B[k] = azB
            each_a1B[k] = bzB
            each_a2B[k] = czB
            each_a3B[k] = dzB
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

            # Accumulating the average and storing the value of the coefficients
            av_a0A += azA
            av_a1A += bzA
            av_a2A += czA
            av_a3A += dzA
            av_a4A += ezA
            av_a0B += azB
            av_a1B += bzB
            av_a2B += czB
            av_a3B += dzB
            av_a4B += ezB
            each_a0A[k] = azA
            each_a1A[k] = bzA
            each_a2A[k] = czA
            each_a3A[k] = dzA
            each_a4A[k] = ezA
            each_a0B[k] = azB
            each_a1B[k] = bzB
            each_a2B[k] = czB
            each_a3B[k] = dzB
            each_a4B[k] = ezB
        
        crossingbeta, crossingU = findintersection(fbetas, fA, fB)
        intersectionbeta_av += crossingbeta
        intersectionbetas[k] = crossingbeta
        temp = converttemp(crossingbeta)
        intersectiontemp_av += temp
        intersectiontemps[k] = temp
        intersectionU_av += crossingU
        intersectionUs[k] = crossingU

        # Finding the exponent this time:
        # Renaming quantities to make eq.s more readable
        x = crossingbeta;
        a1A = each_a1A[k]; a2A = each_a2A[k]; a3A = each_a3A[k]; a4A = each_a4A[k]
        a1B = each_a1B[k]; a2B = each_a2B[k]; a3B = each_a3B[k]; a4B = each_a4B[k]
        fAo = a1A + 2*a2A*x + 3*a3A*x**2 + 4*a4A*x**3
        fBo = a1B + 2*a2B*x + 3*a3B*x**2 + 4*a4B*x**3
        # Setting another quanity I am going to use in both equations
        deltaLdL = float(deltaL)/LA
        log1pdeltaLdL = log(1+deltaLdL) # This must be positive
        # Finding the estimate of the inverse exponent
        inverseexponent = log(abs(fBo/fAo))/log1pdeltaLdL       ####### This will be returned
        exponent1dnu_av += inverseexponent
        exponent1dnus[k] = inverseexponent # Hope this works...
        

        #print "len(fbetas) =", len(fbetas) 
        #print fA
        #print "len(fA) =", len(fA) 
        #print "len(fB) =", len(fB)  
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
            #legend(loc="lower left")
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

    # For A
    av_a0A = av_a0A/no_of_BootstraprunsA
    stddev_a0A = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a0A += (av_a0A-each_a0A[i])*(av_a0A-each_a0A[i])
    stddev_a0A = sqrt(stddev_a0A/(no_of_BootstraprunsA-1))

    av_a1A = av_a1A/no_of_BootstraprunsA
    stddev_a1A = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a1A += (av_a1A-each_a1A[i])*(av_a1A-each_a1A[i])
    stddev_a1A = sqrt(stddev_a1A/(no_of_BootstraprunsA-1))

    av_a2A = av_a2A/no_of_BootstraprunsA
    stddev_a2A = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a2A += (av_a2A-each_a2A[i])*(av_a2A-each_a2A[i])
    stddev_a2A = sqrt(stddev_a2A/(no_of_BootstraprunsA-1))

    av_a3A = av_a3A/no_of_BootstraprunsA
    stddev_a3A = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a3A += (av_a3A-each_a3A[i])*(av_a3A-each_a3A[i])
    stddev_a3A = sqrt(stddev_a3A/(no_of_BootstraprunsA-1))

    av_a4A = av_a4A/no_of_BootstraprunsA
    stddev_a4A = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a4A += (av_a4A-each_a4A[i])*(av_a4A-each_a4A[i])
    stddev_a4A = sqrt(stddev_a4A/(no_of_BootstraprunsA-1))

    # For B
    av_a0B = av_a0B/no_of_BootstraprunsA
    stddev_a0B = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a0B += (av_a0B-each_a0B[i])*(av_a0B-each_a0B[i])
    stddev_a0B = sqrt(stddev_a0B/(no_of_BootstraprunsA-1))

    av_a1B = av_a1B/no_of_BootstraprunsA
    stddev_a1B = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a1B += (av_a1B-each_a1B[i])*(av_a1B-each_a1B[i])
    stddev_a1B = sqrt(stddev_a1B/(no_of_BootstraprunsA-1))

    av_a2B = av_a2B/no_of_BootstraprunsA
    stddev_a2B = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a2B += (av_a2B-each_a2B[i])*(av_a2B-each_a2B[i])
    stddev_a2B = sqrt(stddev_a2B/(no_of_BootstraprunsA-1))

    av_a3B = av_a3B/no_of_BootstraprunsA
    stddev_a3B = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a3B += (av_a3B-each_a3B[i])*(av_a3B-each_a3B[i])
    stddev_a3B = sqrt(stddev_a3B/(no_of_BootstraprunsA-1))

    av_a4B = av_a4B/no_of_BootstraprunsA
    stddev_a4B = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_a4B += (av_a4B-each_a4B[i])*(av_a4B-each_a4B[i])
    stddev_a4B = sqrt(stddev_a4B/(no_of_BootstraprunsA-1))

    coeffsA     = array([av_a0A, av_a1A, av_a2A, av_a3A, av_a4A])
    coeffs_stdA = array([stddev_a0A, stddev_a1A, stddev_a2A, stddev_a3A, stddev_a4A])

    coeffsB     = array([av_a0B, av_a1B, av_a2B, av_a3B, av_a4B])
    coeffs_stdB = array([stddev_a0B, stddev_a1B, stddev_a2B, stddev_a3B, stddev_a4B])

    # Coefficient straight from the Bootstrap method
    exponent1dnu_av = exponent1dnu_av/no_of_BootstraprunsA
    stddev_exponent1dnu = 0
    for i in range(0, no_of_BootstraprunsA):
        stddev_exponent1dnu += (exponent1dnu_av-exponent1dnus[i])*(exponent1dnu_av-exponent1dnus[i])
    stddev_exponent1dnu = sqrt(stddev_exponent1dnu/(no_of_BootstraprunsA-1))

    #exponent1dnus[k]
    # For testing
    '''
    if LA==6:
        print "min(a0) = ", min(each_a0A)
        print "max(a0) = ", max(each_a0A)
        print "min(a1) = ", min(each_a1A)
        print "max(a1) = ", max(each_a1A)
        print "min(a2) = ", min(each_a2A)
        print "max(a2) = ", max(each_a2A)
        print "min(a3) = ", min(each_a3A)
        print "max(a3) = ", max(each_a3A)'''
        
    return intersectionbeta_av, intersectionbeta_stddev, intersectiontemp_av, intersectiontemp_stddev, intersectionU_av, intersectionU_stddev, coeffsA, coeffs_stdA, coeffsB, coeffs_stdB, exponent1dnu_av, stddev_exponent1dnu


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
filenameA = "fcc6x6x6yopen_verylongrunforBinder_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameB = "fcc8x8x8yopen_verylongrunforBinder_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameC = "fcc10x10x10yopen_verylongrunforBinder_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameD = "fcc12x12x12yopen_verylongrunforBinder_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameE = "fcc14x14x14yopen_verylongrunforBinder_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameF = "fcc16x16x16yopen_verylongrunforBinder_Jensen_nnJyz1p04_nnJxy0p3_nnJxzm0p11_nnnJy0p67_nnnJzm0p05_sianDx0p34_Dy1p82_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
'''

# Li
filenameA = "fcc6x6x6yopen_verylongrunforBinder_Li_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameB = "fcc8x8x8yopen_verylongrunforBinder_Li_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameC = "fcc10x10x10yopen_verylongrunforBinder_Li_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameD = "fcc12x12x12yopen_verylongrunforBinder_Li_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameE = "fcc14x14x14yopen_verylongrunforBinder_Li_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"
filenameF = "fcc16x16x16yopen_verylongrunforBinder_Li_nnJyz0p94_nnJxy0p26_nnJxzm0p16_nnnJy0p59_nnnJzm0p11_sianDx0p34_Dy1p92_eq10000_mc1000000_bins100_divseeds_slowcool_binavgs.txt"

##################                                   #######################
bool16 = 1 #0 - we have only up to L=14 / anythong else- we have L=16 too

quadraticnotcubicfit = 2 # 0 gives a line fit to the Binder cumulants, 1 gives a quadratic fit, 2 gives cubic, 3 gives quadruple
toplotornot = 0 # 0, don't plot, anything else leads to a plot

# If we only want to look at a small temperature interval
cutit = 1 # 0: Find a fit for the whole beta range; 1: Find a fit for a small interval
# Jensen
'''
cutlowerAB = 0.778
cutupperAB = 0.786

cutlowerBC = 0.783
cutupperBC = 0.789

cutlowerCD = 0.784
cutupperCD = 0.789

cutlower_rest = 0.784
cutupper_rest = 0.790
'''

# Li
#'''
cutlowerAB = 0.768
cutupperAB = 0.772

cutlowerBC = 0.770
cutupperBC = 0.780

cutlowerCD = 0.772
cutupperCD = 0.776

cutlowerDE = 0.773
cutupperDE = 0.778

cutlower_rest = 0.770
cutupper_rest = 0.780
#'''

#cutlower1 = cutlower2
#cutupper1 = cutupper2

deltaL = 2 # Difference in the size of the systems for the graphs we use for the intersection

print " "
print "L, L+2:"

ib_avA, ib_stdvA, iT_avA, iT_stdvA, iU_avA, iU_stdvA, coeffsA6, coeffs_stdA6, coeffsB8, coeffs_stdB8, exponent1dnu_avA6B8, stddev_exponent1dnuA6B8 = extract_crossing(LA, LB, filenameA, filenameB, quadraticnotcubicfit, toplotornot, cutit, cutlowerAB, cutupperAB, deltaL)
ib_avB, ib_stdvB, iT_avB, iT_stdvB, iU_avB, iU_stdvB, coeffsA8, coeffs_stdA8, coeffsB10, coeffs_stdB10, exponent1dnu_avA8B10, stddev_exponent1dnuA8B10 = extract_crossing(LB, LC, filenameB, filenameC, quadraticnotcubicfit, toplotornot, cutit, cutlowerBC, cutupperBC, deltaL)
ib_avC, ib_stdvC, iT_avC, iT_stdvC, iU_avC, iU_stdvC, coeffsA10, coeffs_stdA10, coeffsB12, coeffs_stdB12, exponent1dnu_avA10B12, stddev_exponent1dnuA10B12 = extract_crossing(LC, LD, filenameC, filenameD, quadraticnotcubicfit, toplotornot, cutit, cutlowerCD, cutupperCD, deltaL)
ib_avD, ib_stdvD, iT_avD, iT_stdvD, iU_avD, iU_stdvD, coeffsA12, coeffs_stdA12, coeffsB14, coeffs_stdB14, exponent1dnu_avA12B14, stddev_exponent1dnuA12B14 = extract_crossing(LD, LE, filenameD, filenameE, quadraticnotcubicfit, toplotornot, cutit, cutlowerDE, cutupperDE, deltaL)
if bool16!=0:
    ib_avE, ib_stdvE, iT_avE, iT_stdvE, iU_avE, iU_stdvE, coeffsA14, coeffs_stdA14, coeffsB16, coeffs_stdB16, exponent1dnu_avA14B16, stddev_exponent1dnuA14B16 = extract_crossing(LE, LF, filenameE, filenameF, quadraticnotcubicfit, toplotornot, cutit, cutlower_rest, cutupper_rest, deltaL)

crossingbetas      = [ib_avA, ib_avB, ib_avC, ib_avD]
crossingbetas_stdv = [ib_stdvA, ib_stdvB, ib_stdvC, ib_stdvD]
crossingtemps      = [iT_avA, iT_avB, iT_avC, iT_avD]
crossingtemps_stdv = [iT_stdvA, iT_stdvB, iT_stdvC, iT_stdvD]
crossingUs         = [iU_avA, iU_avB, iU_avC, iU_avD]
crossingUs_stdv    = [iU_stdvA, iU_stdvB, iU_stdvC, iU_stdvD]
vec1dL             = [1./LA, 1./LB, 1./LC, 1./LD]
vecL               = [LA, LB, LC, LD]
exponent1dnu_array = [exponent1dnu_avA6B8, exponent1dnu_avA8B10, exponent1dnu_avA10B12, exponent1dnu_avA12B14]
exponent1dnu_stdvs = [stddev_exponent1dnuA6B8, stddev_exponent1dnuA8B10, stddev_exponent1dnuA10B12, stddev_exponent1dnuA12B14]

if bool16!=0:
    crossingbetas.append(ib_avE)
    crossingbetas_stdv.append(ib_stdvE)
    crossingtemps.append(iT_avE)
    crossingtemps_stdv.append(iT_stdvE)
    crossingUs.append(iU_avE)
    crossingUs_stdv.append(iU_stdvE)
    vec1dL.append(1./LE)
    vecL.append(LE)
    exponent1dnu_array.append(exponent1dnu_avA14B16)
    exponent1dnu_stdvs.append(stddev_exponent1dnuA14B16)

crossingbetas      = array(crossingbetas)
crossingbetas_stdv = array(crossingbetas_stdv)
crossingtemps      = array(crossingtemps)
crossingtemps_stdv = array(crossingtemps_stdv)
crossingUs         = array(crossingUs)
crossingUs_stdv    = array(crossingUs_stdv)
vec1dL             = array(vec1dL)
vecL               = array(vecL)
exponent1dnu_array = array(exponent1dnu_array)
exponent1dnu_stdvs = array(exponent1dnu_stdvs)

print "Crossing points of graphs"
for betapoint in crossingbetas:
    print betapoint

checkcoeffs = 1
if checkcoeffs==1:
    print "Coeffs of L=6:"
    for i in range(len(coeffsA6)):
        print "coeff[",i,"] = ", coeffsA6[i], "; coeff_std[",i,"] = ", coeffs_stdA6[i]

# Writing to file for experimentation
lent = len(vec1dL)
outfilecrt = open("LandLp2_crtemps_data_Li_001.txt", 'w')
outfilecrU = open("LandLp2_crUs_data_Li_001.txt", 'w')
for i in range(lent):
    outfilecrt.write('%.16f %.16f \n'%(crossingtemps[i], crossingtemps_stdv[i])) # Guess I don't need the standard deviations...
    outfilecrU.write('%.16f %.16f \n'%(crossingUs[i], crossingUs_stdv[i])) # Guess I don't need the standard deviations...
outfilecrt.close()
outfilecrU.close()

maxbeta_inplot = (max(vec1dL)*1.02)  # Misleading name...

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL, crossingbetas, yerr=crossingbetas_stdv, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL, crossingtemps, yerr=crossingtemps_stdv, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'T[K]', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, min(crossingtemps)*0.999, max(crossingtemps)*1.001])
show()


# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL, crossingUs, yerr=crossingUs_stdv, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, (min(crossingUs)-max(crossingUs_stdv))*0.999, (max(crossingUs)+max(crossingUs_stdv))*1.001])
show()

'''
# From regular data points
tempfitparameters, pcov = curve_fit(fittingfunc,  vec1dL,  crossingtemps)#, p0=(guessa,guessb,guessc))
Tc = tempfitparameters[0]; weight = tempfitparameters[1]; exponent = tempfitparameters[2]

print "Parameters of exponential expression: Tc = ", Tc, " weight = ", weight, " exponent = ", exponent

fit_exp = give_fittingfunc(vec1dL, Tc, weight, exponent)
print "Finite size scaling analysis, crossing:", Tc, "in the units of beta:", converttemp(Tc)
print "Weight:", weight
print "Exponent, 1/nu+omega:", exponent
perrstd = sqrt(diag(pcov))
covTalpha    = pcov.item(0,1)
covTexp      = pcov.item(0,2)
covbexpalpha = pcov.item(1,2)
print "Standard deviations, that order:", perrstd[0], ",", perrstd[1], ",", perrstd[2], "."
print "Covariance, ULc and weight:", covTalpha
print "Covariance, ULc and omega:", covTexp
print "Covariance, omega and weight:", covbexpalpha

n = 100
T2weight_min = weight-perrstd[1]
T2weight_max = weight+perrstd[1]
#deltaT2weight = (T2weight_max-T2weight_min)/(n-100)
T2weight_array = linspace(T2weight_min,T2weight_max, n) # Might as well do this
Tcestimates_vs_weigths = zeros(n)
exponentestimates_vs_weigths = zeros(n)
Tcestimates_vs_weigths_stdv = zeros(n)
exponentestimates_vs_weigths_stdv = zeros(n)
for i in range(n):
    weight = T2weight_array[i]
    check_regfit = lambda x,crossing,exponent: crossing + weight*x**(exponent)
    tempfitparameters, pcov = curve_fit(check_regfit, vec1dL, crossingtemps)#, p0=(guessa,guessb,guessc))
    Tc = tempfitparameters[0]; exponent = tempfitparameters[1]
    perrstd = sqrt(diag(pcov))
    Tcestimates_vs_weigths[i] = Tc
    exponentestimates_vs_weigths[i] = exponent
    Tcestimates_vs_weigths_stdv[i] = perrstd[0]
    exponentestimates_vs_weigths_stdv[i] = perrstd[1]

print "Smallest weight:", T2weight_array[0]
print "Uc in range [", (Tcestimates_vs_weigths[0]-Tcestimates_vs_weigths_stdv[0]), ",", (Tcestimates_vs_weigths[0]+Tcestimates_vs_weigths_stdv[0]), "]"
print "1/nu+omega in range [", (exponentestimates_vs_weigths[0]-exponentestimates_vs_weigths_stdv[0]), ",", (exponentestimates_vs_weigths[0]+exponentestimates_vs_weigths_stdv[0]), "]"
print "Largest weight:", T2weight_array[n-1]
print "Tc in range [", (Tcestimates_vs_weigths[n-1]-Tcestimates_vs_weigths_stdv[n-1]), ",", (Tcestimates_vs_weigths[n-1]+Tcestimates_vs_weigths_stdv[n-1]), "]"
print "1/nu+omega in range [", (exponentestimates_vs_weigths[n-1]-exponentestimates_vs_weigths_stdv[n-1]), ",", (exponentestimates_vs_weigths[n-1]+exponentestimates_vs_weigths_stdv[n-1]), "]"

figure(figsize=(6,5))
plot(T2weight_array, Tcestimates_vs_weigths, label='mid')
hold('on')
plot(T2weight_array,Tcestimates_vs_weigths+Tcestimates_vs_weigths_stdv, label='mid+stdv')
plot(T2weight_array, Tcestimates_vs_weigths-Tcestimates_vs_weigths_stdv, label='mid-stdv')
title(r'How $T_c$ varies with the weight', fontsize=14)
xlabel(r'$\alpha$', fontsize=20)
ylabel(r'$T_c$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower right")
show()

figure(figsize=(6,5))
plot(T2weight_array, exponentestimates_vs_weigths, label='mid')
hold('on')
plot(T2weight_array,exponentestimates_vs_weigths+exponentestimates_vs_weigths_stdv, label='mid+stdv')
plot(T2weight_array, exponentestimates_vs_weigths-exponentestimates_vs_weigths_stdv, label='mid-stdv')
title(r'How $\frac{1}{\nu}+\omega$ varies with the weight', fontsize=14)
xlabel(r'$\alpha$', fontsize=20)
ylabel(r'$\frac{1}{\nu}+\omega$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower right")
show()


# To get an error estimate on the critical temperature:
lent = len(crossingtemps)
crossingtemps_max = zeros(lent)
crossingtemps_min = zeros(lent)
for i in range(lent):
    crossingtemps_max[i] = crossingtemps[i]+crossingtemps_stdv[i]
    crossingtemps_min[i] = crossingtemps[i]-crossingtemps_stdv[i]

# From max ends of data points
tempfitparameters, pcov = curve_fit(fittingfunc,  vec1dL,  crossingtemps_max)#, p0=(Tc,weight,exponent))
Tcfvmax = tempfitparameters[0]; wfvmax = tempfitparameters[1]; expfvmax = tempfitparameters[2]
perrstdmax = sqrt(diag(pcov))
covTalphamax    = pcov.item(0,1)
covTexpmax      = pcov.item(0,2)
covbexpalphamax = pcov.item(1,2)
fit_exp_max = give_fittingfunc(vec1dL, Tcfvmax, wfvmax, expfvmax)
print "Adding the stdv:"
print "Finite size scaling analysis, crossing:", Tcfvmax, "in the units of beta:", converttemp(Tcfvmax)
print "Weight:", wfvmax
print "Exponent, 1/nu+omega:", expfvmax
print "Standard deviations, that order:", perrstdmax[0], ",", perrstdmax[1], ",", perrstdmax[2], "."
print "Covariance, ULc and weight:", covTalphamax
print "Covariance, ULc and omega:", covTexpmax
print "Covariance, omega and weight:", covbexpalphamax

# From min ends of data points
tempfitparameters, pcov = curve_fit(fittingfunc,  vec1dL,  crossingtemps_min)#, p0=(Tc,weight,exponent))
Tcfvmin = tempfitparameters[0]; wfvmin = tempfitparameters[1]; expfvmin = tempfitparameters[2]
perrstdmin = sqrt(diag(pcov))
covTalphamin    = pcov.item(0,1)
covTexpmin      = pcov.item(0,2)
covbexpalphamin = pcov.item(1,2)
fit_exp_min = give_fittingfunc(vec1dL, Tcfvmin, wfvmin, expfvmin)
print "Subtracting the stdv:"
print "Finite size scaling analysis, crossing:", Tcfvmin, "in the units of beta:", converttemp(Tcfvmin)
print "Weight:", wfvmin
print "Exponent, 1/nu+omega:", expfvmin
print "Standard deviations, that order:", perrstdmin[0], ",", perrstdmin[1], ",", perrstdmin[2], "."
print "Covariance, ULc and weight:", covTalphamin
print "Covariance, ULc and omega:", covTexpmin
print "Covariance, omega and weight:", covbexpalphamin

figure(figsize=(6,5))
errorbar(vec1dL, crossingtemps, yerr=crossingtemps_stdv, fmt=None, capsize=2)
hold('on')
plot(vec1dL, fit_exp, label='Average fit')
plot(vec1dL, fit_exp_max, label='Max stdv fit')
plot(vec1dL, fit_exp_min, label='Min stdv fit')
title(r'Crossings between Binder cumulant graphs L and L+2', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'T[K]', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="upper center")
show()

# Hardcoding it into an equation
vectbcmtc = crossingtemps-Tc # Datapoints - Tc that I found.
rhs = zeros(len(vectbcmtc))
for i in range(len(rhs)):
    rhs[i] = weight*vec1dL[i]**exponent # what the fit wants me to equate it to
# This should be a line. If its not, I guess the fit was bad.
# log-log plot
figure(figsize=(6,5))
loglog(rhs, vectbcmtc)
hold('on')
loglog(rhs, vectbcmtc, 'o')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'-$\left(\frac{1}{\nu}+\omega\right)$log($L$)+log$\alpha$', fontsize=20)
ylabel(r'log($T^*-T_C$)', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()

print " "

### Fitting U

guessa = 1
guessb = 1
guessc = 0.41
Ufitparameters, pcov = curve_fit(fittingfunc,  vec1dL, crossingUs)#, p0=(guessa,guessb,guessc))
ULc = Ufitparameters[0]; b = Ufitparameters[1]; omega = Ufitparameters[2]
perrUstd = sqrt(diag(pcov))
covUb = pcov.item(0,1)
covUo = pcov.item(0,2)
covbo = pcov.item(1,2)
Ufit_exp = give_fittingfunc(vec1dL, ULc, b, omega)
print "Finite size scaling, average:", ULc
print "Weight:", b
print "Exponent, omega:", omega
print "Standard deviations, that order:", perrUstd[0], ",", perrUstd[1], ",", perrUstd[2], "."
print "Covariance, ULc and weight:", covUb
print "Covariance, ULc and omega:", covUo
print "Covariance, omega and weight:", covbo

n = 100
U2weight_min = b-perrUstd[1]
U2weight_max = b+perrUstd[1]
#deltaT2weight = (T2weight_max-T2weight_min)/(n-100)
U2weight_array = linspace(U2weight_min,U2weight_max, n) # Might as well do this
Ucestimates_vs_weigths = zeros(n)
omegaestimates_vs_weigths = zeros(n)
Ucestimates_vs_weigths_stdv = zeros(n)
omegaestimates_vs_weigths_stdv = zeros(n)
for i in range(n):
    weight = U2weight_array[i]
    check_regfitU = lambda x,crossing,exponent: crossing + weight*x**(exponent)
    Ufitparameters, pcov = curve_fit(check_regfitU, vec1dL, crossingUs)#, p0=(guessa,guessb,guessc))
    Uc = Ufitparameters[0]; exponent = Ufitparameters[1];
    perrUstd = sqrt(diag(pcov))
    Ucestimates_vs_weigths[i] = Uc
    omegaestimates_vs_weigths[i] = exponent
    Ucestimates_vs_weigths_stdv[i] = perrUstd[0]
    omegaestimates_vs_weigths_stdv[i] = perrUstd[1]

print "Smallest weight:", U2weight_array[0]
print "Uc in range [", (Ucestimates_vs_weigths[0]-Ucestimates_vs_weigths_stdv[0]), ",", (Ucestimates_vs_weigths[0]+Ucestimates_vs_weigths_stdv[0]), "]"
print "1/nu+omega in range [", (omegaestimates_vs_weigths[0]-omegaestimates_vs_weigths_stdv[0]), ",", (omegaestimates_vs_weigths[0]+omegaestimates_vs_weigths_stdv[0]), "]"
print "Largest weight:", U2weight_array[n-1]
print "Uc in range [", (Ucestimates_vs_weigths[n-1]-Ucestimates_vs_weigths_stdv[n-1]), ",", (Ucestimates_vs_weigths[n-1]+Ucestimates_vs_weigths_stdv[n-1]), "]"
print "1/nu+omega in range [", (omegaestimates_vs_weigths[n-1]-omegaestimates_vs_weigths_stdv[n-1]), ",", (omegaestimates_vs_weigths[n-1]+omegaestimates_vs_weigths_stdv[n-1]), "]"

figure(figsize=(6,5))
plot(U2weight_array,Ucestimates_vs_weigths+Ucestimates_vs_weigths_stdv, label='mid+stdv')
plot(U2weight_array, Ucestimates_vs_weigths-Ucestimates_vs_weigths_stdv, label='mid-stdv')
hold('on')
plot(U2weight_array, Ucestimates_vs_weigths, label='mid')
title(r'How $U_c$ varies with the weight', fontsize=14)
xlabel(r'$\alpha$', fontsize=20)
ylabel(r'$U_c$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="upper left")
show()

figure(figsize=(6,5))
plot(U2weight_array,omegaestimates_vs_weigths+omegaestimates_vs_weigths_stdv, label='mid+stdv')
plot(U2weight_array, omegaestimates_vs_weigths-omegaestimates_vs_weigths_stdv, label='mid-stdv')
hold('on')
plot(U2weight_array, omegaestimates_vs_weigths, label='mid')
title(r'How $\omega$ varies with the weight', fontsize=14)
xlabel(r'$\alpha$', fontsize=20)
ylabel(r'$\omega$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower center")
show()


# To get an error estimate on the critical temperature:
lent = len(crossingUs)
crossingUs_max = zeros(lent)
crossingUs_min = zeros(lent)
for i in range(lent):
    crossingUs_max[i] = crossingUs[i]+crossingUs_stdv[i]
    crossingUs_min[i] = crossingUs[i]-crossingUs_stdv[i]

Ufitparameters, pcov = curve_fit(fittingfunc,  vec1dL,  crossingUs_max)#, p0=(ULc,b,omega))
ULcmax = Ufitparameters[0]; bmax = Ufitparameters[1]; omegamax = Ufitparameters[2]
perrUstdmax = sqrt(diag(pcov))
Ufit_exp_max = give_fittingfunc(vec1dL, ULcmax, bmax, omegamax)
covUbmax = pcov.item(0,1)
covUomax = pcov.item(0,2)
covbomax = pcov.item(1,2)
print "Adding the standard deviation:"
print "Finite size scaling, average+stdv:", ULcmax
print "Weight:", bmax
print "Exponent, omega:", omegamax
print "Standard deviations, that order:", perrUstdmax[0], ",", perrUstdmax[1], ",", perrUstdmax[2], "."
print "Covariance, ULc and weight:", covUbmax
print "Covariance, ULc and omega:", covUomax
print "Covariance, omega and weight:", covbomax


Ufitparameters, pcov = curve_fit(fittingfunc,  vec1dL,  crossingUs_min)#, p0=(ULc,b,omega))
ULcmin = Ufitparameters[0]; bmin = Ufitparameters[1]; omegamin = Ufitparameters[2]
perrUstdmin = sqrt(diag(pcov))
covUbmin = pcov.item(0,1)
covUomin = pcov.item(0,2)
covbomin = pcov.item(1,2)
Ufit_exp_min = give_fittingfunc(vec1dL, ULcmin, bmin, omegamin)
print "Subtracting the standard deviation:"
print "Finite size scaling, average-stdv:", ULcmin
print "Weight:", bmin
print "Exponent, omega:", omegamin
print "Standard deviations, that order:", perrUstdmin[0], ",", perrUstdmin[1], ",", perrUstdmin[2], "."
print "Covariance, ULc and weight:", covUbmin
print "Covariance, ULc and omega:", covUomin
print "Covariance, omega and weight:", covbomin

figure(figsize=(6,5))
errorbar(vec1dL, crossingUs, yerr=crossingUs_stdv, fmt=None, capsize=2)
hold('on')
plot(vec1dL, Ufit_exp, label='Average')
plot(vec1dL, Ufit_exp_max, label='Max')
plot(vec1dL, Ufit_exp_min, label='Min')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower left")
show()

# Hardcoding it into an equation
vecUbcmUc = abs(crossingUs-ULc) # Datapoints - Tc that I found.
rhsU = zeros(len(vecUbcmUc))
for i in range(len(rhsU)):
    rhsU[i] = abs(b*vec1dL[i]**omega) # what the fit wants me to equate it to
# This should be a line. If its not, I guess the fit was bad.
# log-log plot
figure(figsize=(6,5))
loglog(rhsU, vecUbcmUc)
hold('on')
loglog(rhsU, vecUbcmUc, 'o')
title(r'Crossings between Binder cumulant graphs L and L+2 ', fontsize=14)
xlabel(r'-$\left|\omega\log(L)+\log b\right|$', fontsize=20)
ylabel(r'log($\left|U^*-U_C\right|$)', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()
'''
# Do something like crossingbetas+errorbars, crossingbetas-errorbars etc to get even better errorbars?
'''
# Finding 1/nu from the slope of the cumulant (the derivative of U) at Tc. I have Tc from the curve_fit procedure

# coeffsA6, coeffs_stdA6, coeffsB8, coeffs_stdB8
# Finding the slopes at the crossing
# First pair
slopeA6 = givevalueofderivative(crossingbetas[0], coeffsA6)
slopeB8 = givevalueofderivative(crossingbetas[0], coeffsB8)
# Second pair
slopeA8 = givevalueofderivative(crossingbetas[1], coeffsA8)
slopeB10 = givevalueofderivative(crossingbetas[1], coeffsB10)
# Third pair
slopeA10 = givevalueofderivative(crossingbetas[2], coeffsA10)
slopeB12 = givevalueofderivative(crossingbetas[2], coeffsB12)
# Fourth pair
slopeA12 = givevalueofderivative(crossingbetas[3], coeffsA12)
slopeB14 = givevalueofderivative(crossingbetas[3], coeffsB14)
# Fitfth pair
slopeA14 = givevalueofderivative(crossingbetas[4], coeffsA14)
slopeB16 = givevalueofderivative(crossingbetas[4], coeffsB16)

# Finding the approximations to 1/nu from these data

nuinverse_6_8,   nuinverse_6_8_stdv   = find_1dnu_L(LA, deltaL, crossingbetas[0], crossingbetas_stdv[0], coeffsA6, coeffs_stdA6, coeffsB8, coeffs_stdB8)
nuinverse_8_10,  nuinverse_8_10_stdv  = find_1dnu_L(LB, deltaL, crossingbetas[1], crossingbetas_stdv[1], coeffsA8, coeffs_stdA8, coeffsB10, coeffs_stdB10)
nuinverse_10_12, nuinverse_10_12_stdv = find_1dnu_L(LC, deltaL, crossingbetas[2], crossingbetas_stdv[2], coeffsA10, coeffs_stdA10, coeffsB12, coeffs_stdB12)
nuinverse_12_14, nuinverse_12_14_stdv = find_1dnu_L(LD, deltaL, crossingbetas[3], crossingbetas_stdv[3], coeffsA12, coeffs_stdA12, coeffsB14, coeffs_stdB14)
nuinverse_14_16, nuinverse_14_16_stdv = find_1dnu_L(LE, deltaL, crossingbetas[4], crossingbetas_stdv[4], coeffsA14, coeffs_stdA14, coeffsB16, coeffs_stdB16)

nuinverse_array = array([nuinverse_6_8, nuinverse_8_10, nuinverse_10_12, nuinverse_12_14, nuinverse_14_16])
nuinverse_stdv_array = array([nuinverse_6_8_stdv, nuinverse_8_10_stdv, nuinverse_10_12_stdv, nuinverse_12_14_stdv, nuinverse_14_16_stdv])

# Plotting
figure(figsize=(6,5))
plot(vecL, nuinverse_array, 'o')
hold('on')
#errorbar(vecL, nuinverse_array, nuinverse_stdv_array, fmt=None, capsize=2)
xlabel('L', fontsize=20)
ylabel(r'1/$\nu$', fontsize=20)
title(r'Estimation of 1/$\nu$, the dumb way', fontsize=14)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
show()

print "1/nu, from L=14 and L = 16:", nuinverse_14_16, "; standard deviation: ", nuinverse_14_16_stdv


print "Directly from the Bootstrap procedure:"
smallest_points = exponent1dnu_array-exponent1dnu_stdvs
largest_points = exponent1dnu_array+exponent1dnu_stdvs

figure(figsize=(6,5))
errorbar(vecL, exponent1dnu_array, exponent1dnu_stdvs, fmt=None, capsize=2)
hold('on')
plot(vecL, exponent1dnu_array, 'o')
xlabel('L', fontsize=20)
ylabel(r'1/$\nu$', fontsize=20)
title(r'Estimation of 1/$\nu$, L, L+2', fontsize=14)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([min(vecL)-0.1-0.1, max(vecL)+0.1, min(smallest_points)-1, max(largest_points)+1])
show()

leam1 = len(exponent1dnu_array)-1
exponent1dnu_av1416 = exponent1dnu_array[leam1]
stddev_exponent1dnu1416 = exponent1dnu_stdvs[leam1]
print "From L=14 and L=16: 1/nu = ",  exponent1dnu_av1416, "; standard deviation = ", stddev_exponent1dnu1416
'''
#################### BOOTSTRAP L, L+4######################################

quadraticnotcubicfit = 3 # 0 gives a line fit to the Binder cumulants, 2 gives a quadratic fit, 3 gives cubic
toplotornot = 1 # 0, don't plot, anything else leads to a plot
# A-6, B-8, C-10, D-12, E-14, F-16
# 6-10-14, A-C-E
# 8-12-16, B-D-F
# Jensen
'''
cutlowerAC = 0.780
cutupperAC = 0.786

cutlowerBD = 0.784
cutupperBD = 0.789

cutlowerCE = 0.784
cutupperCE = 0.789

cutlower_rest = 0.784
cutupper_rest = 0.789
'''
# Li
#'''
cutlowerAC = 0.770
cutupperAC = 0.775

cutlowerBD = 0.772
cutupperBD = 0.780

cutlowerCE = 0.772
cutupperCE = 0.780

cutlower_rest = 0.770
cutupper_rest = 0.780
#'''


deltaL = 4 # Difference in the size of the systems for the graphs we use for the intersection

print " "
print "L, L+4:"

ib_avA4, ib_stdvA4, iT_avA4, iT_stdvA4, iU_avA4, iU_stdvA4, coeffsA6_L4, coeffs_stdA6_L4, coeffsB10_L4, coeffs_stdB10_L4, exponent1dnu_avA6B10, stddev_exponent1dnuA6B10 = extract_crossing(LA, LC, filenameA, filenameC, quadraticnotcubicfit, toplotornot, cutit, cutlowerAC, cutupperAC, deltaL)
ib_avB4, ib_stdvB4, iT_avB4, iT_stdvB4, iU_avB4, iU_stdvB4, coeffsA8_L4, coeffs_stdA8_L4, coeffsB12_L4, coeffs_stdB12_L4, exponent1dnu_avA8B12, stddev_exponent1dnuA8B12  = extract_crossing(LB, LD, filenameB, filenameD, quadraticnotcubicfit, toplotornot, cutit, cutlowerBD, cutupperBD, deltaL)
ib_avC4, ib_stdvC4, iT_avC4, iT_stdvC4, iU_avC4, iU_stdvC4, coeffsA10_L4, coeffs_stdA10_L4, coeffsB14_L4, coeffs_stdB14_L4, exponent1dnu_avA10B14, stddev_exponent1dnuA10B14  = extract_crossing(LC, LE, filenameC, filenameE, quadraticnotcubicfit, toplotornot, cutit, cutlowerCE, cutupperCE, deltaL)

if bool16!=0:
    ib_avD4, ib_stdvD4, iT_avD4, iT_stdvD4, iU_avD4, iU_stdvD4, coeffsA12_L4, coeffs_stdA12_L4, coeffsB16_L4, coeffs_stdB16_L4, exponent1dnu_avA12B16, stddev_exponent1dnuA12B16  = extract_crossing(LD, LF, filenameD, filenameF, quadraticnotcubicfit, toplotornot, cutit, cutlower_rest, cutupper_rest, deltaL)

crossingbetas4      = [ib_avA4, ib_avB4, ib_avC4]
crossingbetas_stdv4 = [ib_stdvA4, ib_stdvB4, ib_stdvC4]
crossingtemps4      = [iT_avA4, iT_avB4, iT_avC4]
crossingtemps_stdv4 = [iT_stdvA4, iT_stdvB4, iT_stdvC4]
crossingUs4         = [iU_avA4, iU_avB4, iU_avC4]
crossingUs_stdv4    = [iU_stdvA4, iU_stdvB4, iU_stdvC4]
vec1dL4             = [1./LA, 1./LB, 1./LC]
vecLsp4             = [LA, LB, LC]
exponent1dnu_array4 = [exponent1dnu_avA6B10, exponent1dnu_avA8B12, exponent1dnu_avA10B14]
exponent1dnu_stdvs4 = [stddev_exponent1dnuA6B10, stddev_exponent1dnuA8B12, stddev_exponent1dnuA10B14]

if bool16!=0:
    crossingbetas4.append(ib_avD4)
    crossingbetas_stdv4.append(ib_stdvD4)
    crossingtemps4.append(iT_avD4)
    crossingtemps_stdv4.append(iT_stdvD4)
    crossingUs4.append(iU_avD4)
    crossingUs_stdv4.append(iU_stdvD4)
    vec1dL4.append(1./LD)
    vecLsp4.append(LD)
    exponent1dnu_array4.append(exponent1dnu_avA12B16)
    exponent1dnu_stdvs4.append(stddev_exponent1dnuA12B16)


crossingbetas4      = array(crossingbetas4)
crossingbetas_stdv4 = array(crossingbetas_stdv4)
crossingtemps4      = array(crossingtemps4)
crossingtemps_stdv4 = array(crossingtemps_stdv4)
crossingUs4         = array(crossingUs4)
crossingUs_stdv4    = array(crossingUs_stdv4)
vec1dL4             = array(vec1dL4)
vecLsp4             = array(vecLsp4)
exponent1dnu_array4 = array(exponent1dnu_array4)
exponent1dnu_stdvs4 = array(exponent1dnu_stdvs4)

print "Crossing points of graphs (beta)"
for betapoint in crossingbetas4:
    print betapoint

maxbeta_inplot = (max(vec1dL4)*1.02)

# Writing to file for experimentation
lent = len(vec1dL4)
outfilecrt4 = open("LandLp4_crtemps_data_Li_001.txt", 'w')
outfilecrU4 = open("LandLp4_crUs_data_Li_001.txt", 'w')
for i in range(lent):
    outfilecrt4.write('%.16f %.16f \n'%(crossingtemps4[i], crossingtemps_stdv4[i])) # Guess I don't need the standard deviations...
    outfilecrU4.write('%.16f %.16f \n'%(crossingUs4[i], crossingUs_stdv4[i])) # Guess I don't need the standard deviations...
outfilecrt4.close()
outfilecrU4.close()

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL4, crossingbetas4, yerr=crossingbetas_stdv4, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+4 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$\beta$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, min(crossingbetas4)*0.999, max(crossingbetas4)*1.001])
show()

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL4, crossingtemps4, yerr=crossingtemps_stdv4, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+4 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'T[K]', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, min(crossingtemps4)*0.999, max(crossingtemps4)*1.001])
show()

# Regular plot
figure(figsize=(6,5))
errorbar(vec1dL4, crossingUs4, yerr=crossingUs_stdv4, fmt=None, capsize=2)
hold('on')
title(r'Crossings between Binder cumulant graphs L and L+4 ', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([0, maxbeta_inplot, (min(crossingUs4)-max(crossingUs_stdv4))*0.999, (max(crossingUs4)+max(crossingUs_stdv4))*1.001])
show()

'''
guessa = 1
guessb = 1
guessc = 14.5
tempfitparameters, pcov = curve_fit(fittingfunc,  vec1dL4,  crossingtemps4, method='trf')#, p0=(guessa,guessb,guessc))
Tc4 = tempfitparameters[0]; weight4 = tempfitparameters[1]; exponent4 = tempfitparameters[2]
perrstd4 = sqrt(diag(pcov))
print "Parameters of exponential expression: Tc = ", Tc4, " weight = ", weight4, " exponent = ", exponent4
print "Weight:", weight4
print "Exponent, 1/nu+omega:", exponent4
perrstd4 = sqrt(diag(pcov))
print "Standard deviations, that order:", perrstd4[0], ",", perrstd4[1], ",", perrstd4[2], "."
covTalpha4    = pcov.item(0,1)
covTexp4      = pcov.item(0,2)
covbexpalpha4 = pcov.item(1,2)
print "Covariance, ULc and weight:", covTalpha4
print "Covariance, ULc and omega:", covTexp4
print "Covariance, omega and weight:", covbexpalpha4

fit_exp = give_fittingfunc(vec1dL4, Tc4, weight4, exponent4)

n = 100
T4weight_min = weight4-perrstd4[1]
T4weight_max = weight4+perrstd4[1]
#deltaT2weight = (T2weight_max-T2weight_min)/(n-100)
T4weight_array = linspace(T4weight_min,T4weight_max, n) # Might as well do this
Tcestimates_vs_weigths = zeros(n)
exponentestimates_vs_weigths = zeros(n)
Tcestimates_vs_weigths_stdv = zeros(n)
exponentestimates_vs_weigths_stdv = zeros(n)
for i in range(n):
    weight = T2weight_array[i]
    check_regfit = lambda x,crossing,exponent: crossing + weight*x**(exponent)
    tempfitparameters, pcov = curve_fit(check_regfit, vec1dL4, crossingtemps4)#, p0=(guessa,guessb,guessc))
    Tc = tempfitparameters[0]; exponent = tempfitparameters[1]
    perrstd = sqrt(diag(pcov))
    Tcestimates_vs_weigths[i] = Tc
    exponentestimates_vs_weigths[i] = exponent
    Tcestimates_vs_weigths_stdv[i] = perrstd[0]
    exponentestimates_vs_weigths_stdv[i] = perrstd[1]

print "Smallest weight:", T4weight_array[0]
print "Tc in range [", (Tcestimates_vs_weigths[0]-Tcestimates_vs_weigths_stdv[0]), ",", (Tcestimates_vs_weigths[0]+Tcestimates_vs_weigths_stdv[0]), "]"
print "1/nu+omega in range [", (exponentestimates_vs_weigths[0]-exponentestimates_vs_weigths_stdv[0]), ",", (exponentestimates_vs_weigths[0]+exponentestimates_vs_weigths_stdv[0]), "]"
print "Largest weight:", T4weight_array[n-1]
print "Tc in range [", (Tcestimates_vs_weigths[n-1]-Tcestimates_vs_weigths_stdv[n-1]), ",", (Tcestimates_vs_weigths[n-1]+Tcestimates_vs_weigths_stdv[n-1]), "]"
print "1/nu+omega in range [", (exponentestimates_vs_weigths[n-1]-exponentestimates_vs_weigths_stdv[n-1]), ",", (exponentestimates_vs_weigths[n-1]+exponentestimates_vs_weigths_stdv[n-1]), "]"

figure(figsize=(6,5))
plot(T4weight_array,Tcestimates_vs_weigths+Tcestimates_vs_weigths_stdv, label='mid+stdv')
plot(T4weight_array, Tcestimates_vs_weigths-Tcestimates_vs_weigths_stdv, label='mid-stdv')
hold('on')
plot(T4weight_array, Tcestimates_vs_weigths, label='mid')
title(r'How $T_c$ varies with the weight', fontsize=14)
xlabel(r'$\alpha$', fontsize=20)
ylabel(r'$T_c$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower right")
show()

figure(figsize=(6,5))
plot(T4weight_array,exponentestimates_vs_weigths+exponentestimates_vs_weigths_stdv, label='mid+stdv')
plot(T4weight_array, exponentestimates_vs_weigths-exponentestimates_vs_weigths_stdv, label='mid-stdv')
hold('on')
plot(T4weight_array, exponentestimates_vs_weigths, label='mid')
title(r'How $\frac{1}{\nu}+\omega$ varies with the weight', fontsize=14)
xlabel(r'$\alpha$', fontsize=20)
ylabel(r'$\frac{1}{\nu}+\omega$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower right")
show()


# To get an error estimate on the critical temperature:
lent = len(crossingtemps4)
crossingtemps_max4 = zeros(lent)
crossingtemps_min4 = zeros(lent)
for i in range(lent):
    crossingtemps_max4[i] = crossingtemps4[i]+crossingtemps_stdv4[i]
    crossingtemps_min4[i] = crossingtemps4[i]-crossingtemps_stdv4[i]

tempfitparameters, pcov = curve_fit(fittingfunc,  vec1dL4,  crossingtemps_max4, method='trf')#, p0=(Tc4,weight4,exponent4))
Tcmax4 = tempfitparameters[0]; weightmax4 = tempfitparameters[1]; exponentmax4 = tempfitparameters[2]
perrstd4max = sqrt(diag(pcov))
fit_exp_max4 = give_fittingfunc(vec1dL4, Tcmax4, weightmax4, exponentmax4)
print "Adding the stdv:"
print "Finite size scaling analysis, crossing:", Tcmax4, "in the units of beta:", converttemp(Tcmax4)
print "Weight:", weightmax4
print "Exponent, 1/nu+omega:", exponentmax4
print "Standard deviations, that order:", perrstd4max[0], ",", perrstd4max[1], ",", perrstd4max[2], "."
covTalpha4max    = pcov.item(0,1)
covTexp4max      = pcov.item(0,2)
covbexpalpha4max = pcov.item(1,2)
print "Covariance, ULc and weight:", covTalpha4max
print "Covariance, ULc and omega:", covTexp4max
print "Covariance, omega and weight:", covbexpalpha4max


tempfitparameters, pcov = curve_fit(fittingfunc,  vec1dL4,  crossingtemps_min4, method='trf')#, p0=(Tc4,weight4,exponent4))
Tcmin4 = tempfitparameters[0]; weightmin4 = tempfitparameters[1]; exponentmin4 = tempfitparameters[2]
perrstd4min = sqrt(diag(pcov))
fit_exp_min4 = give_fittingfunc(vec1dL4, Tcmin4, weightmin4, exponentmin4)
print "Subtracting the stdv:"
print "Finite size scaling analysis, crossing:", Tcmin4, "in the units of beta:", converttemp(Tcmin4)
print "Weight:", weightmin4
print "Exponent, 1/nu+omega:", exponentmin4
print "Standard deviations, that order:", perrstd4min[0], ",", perrstd4min[1], ",", perrstd4min[2], "."
covTalpha4min    = pcov.item(0,1)
covTexp4min      = pcov.item(0,2)
covbexpalpha4min = pcov.item(1,2)
print "Covariance, ULc and weight:", covTalpha4min
print "Covariance, ULc and omega:", covTexp4min
print "Covariance, omega and weight:", covbexpalpha4min

figure(figsize=(6,5))
errorbar(vec1dL4, crossingtemps4, yerr=crossingtemps_stdv4, fmt=None, capsize=2)
hold('on')
plot(vec1dL4, fit_exp, label='Average fit')
plot(vec1dL4, fit_exp_max4, label='Max stdv fit')
plot(vec1dL4, fit_exp_min4, label='Min stdv fit')
title(r'Crossings between Binder cumulant graphs L and L+4', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'T[K]', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="upper center")
show()

# Hardcoding it into an equation
vectbcmtc4 = crossingtemps4-Tc4 # Datapoints - Tc that I found.
rhs4 = zeros(len(vectbcmtc4))
for i in range(len(rhs4)):
    rhs4[i] = weight4*vec1dL4[i]**exponent4 # what the fit wants me to equate it to
# This should be a line. If its not, I guess the fit was bad.
# log-log plot
figure(figsize=(6,5))
loglog(rhs4, vectbcmtc4)
hold('on')
loglog(rhs4, vectbcmtc4, 'o')
title(r'Crossings between Binder cumulant graphs L and L+4 ', fontsize=14)
xlabel(r'-$\left(\frac{1}{\nu}+\omega\right)$log($L$)+log$\alpha$', fontsize=20)
ylabel(r'log($T^*-T_C$)', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()

print " "

### Fitting U
guessa = 1
guessb = 1
guessc = 0.41
Ufitparameters, pcov = curve_fit(fittingfunc,  vec1dL4,  crossingUs4, method='trf')#, p0=(guessa,guessb,guessc))
ULc4 = Ufitparameters[0]; b4 = Ufitparameters[1]; Uexponent4 = Ufitparameters[2]
perrUstd4 = sqrt(diag(pcov))
covUb = pcov.item(0,1)
covUo = pcov.item(0,2)
covbo = pcov.item(1,2)
Ufit_exp = give_fittingfunc(vec1dL4, ULc4, b4, Uexponent4)
print "Estimate of UL, L,L+4:", ULc4
print "Weight:", b4
print "Exponent, omega:", Uexponent4
print "Standard deviations, that order:", perrUstd4[0], ",", perrUstd4[1], ",", perrUstd4[2], "."
print "Covariance, ULc and weight:", covUb
print "Covariance, ULc and omega:", covUo
print "Covariance, omega and weight:", covbo

n = 100
U4weight_min = b4-perrUstd4[1]
U4weight_max = b4+perrUstd4[1]
#deltaT2weight = (T2weight_max-T2weight_min)/(n-100)
U4weight_array = linspace(U4weight_min,U4weight_max, n) # Might as well do this
Ucestimates_vs_weigths = zeros(n)
omegaestimates_vs_weigths = zeros(n)
Ucestimates_vs_weigths_stdv = zeros(n)
omegaestimates_vs_weigths_stdv = zeros(n)
for i in range(n):
    weight = U4weight_array[i]
    check_regfit = lambda x,crossing,exponent: crossing + weight*x**(exponent)
    Ufitparameters, pcov = curve_fit(check_regfit, vec1dL4, crossingUs4)#, p0=(guessa,guessb,guessc))
    Uc = Ufitparameters[0]; exponent = Ufitparameters[1]
    perrUstd = sqrt(diag(pcov))
    Ucestimates_vs_weigths[i] = Uc
    omegaestimates_vs_weigths[i] = exponent
    Ucestimates_vs_weigths_stdv[i] = perrUstd[0]
    omegaestimates_vs_weigths_stdv[i] = perrUstd[1]

print "Smallest weight:", U4weight_array[0]
print "Uc in range [", (Ucestimates_vs_weigths[0]-Ucestimates_vs_weigths_stdv[0]), ",", (Ucestimates_vs_weigths[0]+Ucestimates_vs_weigths_stdv[0]), "]"
print "1/nu+omega in range [", (omegaestimates_vs_weigths[0]-omegaestimates_vs_weigths_stdv[0]), ",", (omegaestimates_vs_weigths[0]+omegaestimates_vs_weigths_stdv[0]), "]"
print "Largest weight:", U4weight_array[n-1]
print "Tc in range [", (Ucestimates_vs_weigths[n-1]-Ucestimates_vs_weigths_stdv[n-1]), ",", (Ucestimates_vs_weigths[n-1]+Ucestimates_vs_weigths_stdv[n-1]), "]"
print "1/nu+omega in range [", (omegaestimates_vs_weigths[n-1]-omegaestimates_vs_weigths_stdv[n-1]), ",", (omegaestimates_vs_weigths[n-1]+omegaestimates_vs_weigths_stdv[n-1]), "]"

figure(figsize=(6,5))
plot(U4weight_array, Ucestimates_vs_weigths, label='mid')
hold('on')
plot(U4weight_array, Ucestimates_vs_weigths+Ucestimates_vs_weigths_stdv, label='mid+stdv')
plot(U4weight_array, Ucestimates_vs_weigths-Ucestimates_vs_weigths_stdv, label='mid-stdv')
title(r'How $U_c$ varies with the weight', fontsize=14)
xlabel(r'$\alpha$', fontsize=20)
ylabel(r'$U_c$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="upper left")
show()

figure(figsize=(6,5))
plot(U4weight_array, omegaestimates_vs_weigths, label='mid')
hold('on')
plot(U4weight_array,omegaestimates_vs_weigths+omegaestimates_vs_weigths_stdv, label='mid+stdv')
plot(U4weight_array, omegaestimates_vs_weigths-omegaestimates_vs_weigths_stdv, label='mid-stdv')
title(r'How $\omega$ varies with the weight', fontsize=14)
xlabel(r'$\alpha$', fontsize=20)
ylabel(r'$\omega$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower center")
show()

# To get an error estimate on the critical temperature:
lent = len(crossingUs4)
crossingUs_max = zeros(lent)
crossingUs_min = zeros(lent)
for i in range(lent):
    crossingUs_max[i] = crossingUs4[i]+crossingUs_stdv4[i]
    crossingUs_min[i] = crossingUs4[i]-crossingUs_stdv4[i]

Ufitparameters, pcov = curve_fit(fittingfunc,  vec1dL4,  crossingUs_max, method='trf')#, p0=(ULc4,b4,Uexponent4))
ULcmax4 = Ufitparameters[0]; bmax4 = Ufitparameters[1]; Uexponentmax4 = Ufitparameters[2]
perrUstd4max = sqrt(diag(pcov))
covUbmax = pcov.item(0,1)
covUomax = pcov.item(0,2)
covbomax = pcov.item(1,2)
Ufit_exp_max = give_fittingfunc(vec1dL4, ULcmax4, bmax4, Uexponentmax4)
print "Adding the stdv"
print "Estimate of UL, L,L+4:", ULcmax4
print "Weight:", bmax4
print "Exponent, omega:", Uexponentmax4
print "Standard deviations, that order:", perrUstd4max[0], ",", perrUstd4max[1], ",", perrUstd4max[2], "."
print "Covariance, ULc and weight:", covUbmax
print "Covariance, ULc and omega:", covUomax
print "Covariance, omega and weight:", covbomax

Ufitparameters, pcov = curve_fit(fittingfunc,  vec1dL4,  crossingUs_min, method='trf')#, p0=(ULc4,b4,Uexponent4))
ULcmin4 = Ufitparameters[0]; bmin4 = Ufitparameters[1]; Uexponentmin4 = Ufitparameters[2]
perrUstd4min = sqrt(diag(pcov))
covUbmin = pcov.item(0,1)
covUomin = pcov.item(0,2)
covbomin = pcov.item(1,2)
Ufit_exp_min = give_fittingfunc(vec1dL4, ULcmin4, bmin4, Uexponentmin4)
print "Subtracting the stdv"
print "Estimate of UL, L,L+4:", ULcmin4
print "Weight:", bmin4
print "Exponent, omega:", Uexponentmin4
print "Standard deviations, that order:", perrUstd4min[0], ",", perrUstd4min[1], ",", perrUstd4min[2], "."
print "Covariance, ULc and weight:", covUbmin
print "Covariance, ULc and omega:", covUomin
print "Covariance, omega and weight:", covbomin

figure(figsize=(6,5))
errorbar(vec1dL4, crossingUs4, yerr=crossingUs_stdv4, fmt=None, capsize=2)
hold('on')
plot(vec1dL4, Ufit_exp, label='Average')
plot(vec1dL4, Ufit_exp_max, label='Max')
plot(vec1dL4, Ufit_exp_min, label='Min')
title(r'Crossings between Binder cumulant graphs L and L+4', fontsize=14)
xlabel(r'$1/L$', fontsize=20)
ylabel(r'$U_{L,z}$', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
legend(loc="lower left")
show()
'''
# Writing to file
r'''
outfilecrt4 = open("LandLp4_crtemps_data.txt", 'w')
outfilecrU4 = open("LandLp4_crUs_data.txt", 'w')
outfileparameterst4 = open("LandLp4_fitparametersT.txt", 'w')
outfileparametersU4 = open("LandLp4_fitparametersU.txt", 'w')
for i in range(lent):
    outfilecrt4.write('%.16f %.16f \n'%(crossingtemps4[i],crossingtemps_stdv4[i])) # Guess I don't need the standard deviations...
    outfilecrU4.write('%.16f %.16f \n'%(crossingUs4[i],crossingUs_stdv4[i])) # Guess I don't need the standard deviations...
outfileparameterst4.write("Parameters: \n")
outfileparameterst4.write('%.16f '%Tc4)
outfileparameterst4.write('%.16f '%weight4)
outfileparameterst4.write('%.16f \n \n'%exponent4)
outfileparameterst4.write("Standard deviations: \n")
outfileparameterst4.write('%.16f '%perrstd[0])
outfileparameterst4.write('%.16f '%perrstd[1])
outfileparameterst4.write('%.16f '%perrstd[2])

outfileparametersU4.write("Parameters: \n")
outfileparametersU4.write('%.16f '%ULc4)
outfileparametersU4.write('%.16f '%b4)
outfileparametersU4.write('%.16f \n \n'%Uexponent4)
outfileparametersU4.write("Standard deviations: \n")
outfileparametersU4.write('%.16f '%perrUstd4[0])
outfileparametersU4.write('%.16f '%perrUstd4[1])
outfileparametersU4.write('%.16f '%perrUstd4[2])
outfilecrt4.close()
outfilecrU4.close()
outfileparameterst4.close()
outfileparametersU4.close()
'''
'''
# Hardcoding it into an equation
vecUbcmUc4 = abs(crossingUs4-ULc4) # Datapoints - Tc that I found.
rhsU4 = zeros(len(vecUbcmUc4))
for i in range(len(rhsU4)):
    rhsU4[i] = abs(b4*vec1dL4[i]**Uexponent4) # what the fit wants me to equate it to
# This should be a line. If its not, I guess the fit was bad.
# log-log plot
figure(figsize=(6,5))
loglog(rhsU4, vecUbcmUc4)
hold('on')
loglog(rhsU4, vecUbcmUc4, 'o')
title(r'Crossings between Binder cumulant graphs L and L+4 ', fontsize=14)
xlabel(r'$-\left|\omega\log(L)+\log b\right|$', fontsize=20)
ylabel(r'log($\left|U^*-U_C\right|$)', fontsize=20)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
#axis([0, maxbeta_inplot, min(crossingbetas)*0.999, max(crossingbetas)*1.001])
show()
'''
'''
# Finding 1/nu from the slope of the cumulant (the derivative of U) at Tc. I have Tc from the curve_fit procedure
# coeffsA6, coeffs_stdA6, coeffsB8, coeffs_stdB8
# Finding the slopes at the crossing
# First pair
slopeA6_L4  = givevalueofderivative(crossingbetas4[0], coeffsA6_L4)
slopeB10_L4 = givevalueofderivative(crossingbetas4[0], coeffsB10_L4)
# Second pair
slopeA8_L4  = givevalueofderivative(crossingbetas4[1], coeffsA8_L4)
slopeB12_L4 = givevalueofderivative(crossingbetas4[1], coeffsB12_L4)
# Third pair
slopeA10_L4 = givevalueofderivative(crossingbetas4[2], coeffsA10_L4)
slopeB14_L4 = givevalueofderivative(crossingbetas4[2], coeffsB14_L4)
# Fourth pair
slopeA12_L4 = givevalueofderivative(crossingbetas4[3], coeffsA12_L4)
slopeB16_L4 = givevalueofderivative(crossingbetas4[3], coeffsB16_L4)

# Finding the approximations to 1/nu from these data

nuinverse_6_10,   nuinverse_6_10_stdv  = find_1dnu_L(LA, deltaL, crossingbetas4[0], crossingbetas_stdv4[0], coeffsA6_L4, coeffs_stdA6_L4, coeffsB10_L4, coeffs_stdB10_L4)
nuinverse_8_12,   nuinverse_8_12_stdv  = find_1dnu_L(LB, deltaL, crossingbetas4[1], crossingbetas_stdv4[1], coeffsA8_L4, coeffs_stdA8_L4, coeffsB12_L4, coeffs_stdB12_L4)
nuinverse_10_14,  nuinverse_10_14_stdv = find_1dnu_L(LC, deltaL, crossingbetas4[2], crossingbetas_stdv4[2], coeffsA10_L4, coeffs_stdA10_L4, coeffsB14_L4, coeffs_stdB14_L4)
nuinverse_12_16,  nuinverse_12_16_stdv = find_1dnu_L(LD, deltaL, crossingbetas4[3], crossingbetas_stdv4[3], coeffsA12_L4, coeffs_stdA12_L4, coeffsB16_L4, coeffs_stdB16_L4)

nuinverse_array_L4 = array([nuinverse_6_10, nuinverse_8_12, nuinverse_10_14, nuinverse_12_16])
nuinverse_stdv_array_L4 = array([nuinverse_6_10_stdv, nuinverse_8_12_stdv, nuinverse_10_14_stdv, nuinverse_12_16_stdv])

# Plotting
figure(figsize=(6,5))
plot(vecLsp4, nuinverse_array_L4, 'o')
hold('on')
#errorbar(vecLsp4, nuinverse_array_L4, nuinverse_stdv_array_L4, fmt=None, capsize=2)
xlabel('L', fontsize=20)
ylabel(r'1/$\nu$', fontsize=20)
title(r'Estimation of 1/$\nu$, L, L+4, the dumb way', fontsize=14)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
show()

print "1/nu, from L=12 and L = 16:", nuinverse_12_16, "; standard deviation: ", nuinverse_12_16_stdv

print "Directly from the Bootstrap procedure:"
smallest_points = exponent1dnu_array4-exponent1dnu_stdvs4
largest_points = exponent1dnu_array4+exponent1dnu_stdvs4

figure(figsize=(6,5))
errorbar(vecLsp4, exponent1dnu_array4, exponent1dnu_stdvs4, fmt=None, capsize=2)
hold('on')
plot(vecLsp4, exponent1dnu_array4, 'o')
xlabel('L', fontsize=20)
ylabel(r'1/$\nu$', fontsize=20)
title(r'Estimation of 1/$\nu$, L, L+4', fontsize=14)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([min(vecLsp4)-0.1, max(vecLsp4)+0.1, min(smallest_points)-1, max(largest_points)+1])
show()

leam1 = len(exponent1dnu_array4)-1
exponent1dnu_av1216 = exponent1dnu_array4[leam1]
stddev_exponent1dnu1216 = exponent1dnu_stdvs4[leam1]
print "From L=12 and L=16: 1/nu = ", exponent1dnu_av1216 , "; standard deviation = ", stddev_exponent1dnu1216

cutlowerAB = 0.784
cutupperAB = 0.789

cutit = 1 # cutit = 1, cut it, cutit = 0, do nothing.
ib_av1016, ib_stdv1016, iT_av1016, iT_stdv1016, iU_av1016, iU_stdv1016, coeffsA10, coeffs_stdA10, coeffsB16, coeffs_stdB16, exponent1dnu_av1016, stddev_exponent1dnu1016 = extract_crossing(10, 16, filenameC, filenameF, quadraticnotcubicfit, toplotornot, cutit, cutlowerAB, cutupperAB, 6)

print "From L=10 and L=16: 1/nu = ",  exponent1dnu_av1016, "; standard deviation = ", stddev_exponent1dnu1016

ib_av816, ib_stdv816, iT_av816, iT_stdv816, iU_av816, iU_stdv816, coeffsA10, coeffs_stdA10, coeffsB16, coeffs_stdB16, exponent1dnu_av816, stddev_exponent1dnu816 = extract_crossing(8, 16, filenameB, filenameF, quadraticnotcubicfit, toplotornot, cutit, cutlowerAB, cutupperAB, 8)

print "From L=8 and L=16: 1/nu = ",  exponent1dnu_av816, "; standard deviation = ", stddev_exponent1dnu816

deltaL = array([2,4,6,8])
exp_from_deltaLto16       = array([exponent1dnu_av1416,exponent1dnu_av1216,exponent1dnu_av1016,exponent1dnu_av816])
exp_from_deltaLto16_stdvs = array([stddev_exponent1dnu1416,stddev_exponent1dnu1216,stddev_exponent1dnu1016,stddev_exponent1dnu816])

# For plotting:
smallest_points = exp_from_deltaLto16-exp_from_deltaLto16_stdvs
largest_points = exp_from_deltaLto16+exp_from_deltaLto16_stdvs

figure(figsize=(6,5))
errorbar(deltaL, exp_from_deltaLto16, exp_from_deltaLto16_stdvs, fmt=None, capsize=2)
hold('on')
plot(deltaL, exp_from_deltaLto16, 'o')
xlabel('$\Delta L$', fontsize=20)
ylabel(r'1/$\nu$', fontsize=20)
title(r'1/$\nu$  from  $L_1,L_2$ = 16 - $\Delta L$,16 ', fontsize=14)
tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
axis([min(deltaL)-0.1, max(deltaL)+0.1, min(smallest_points)-1, max(largest_points)+1])
show()
print " "
print "deltaL       ", "1/nu       ", "stdv, 1/nu"
for i in range(len(deltaL)):
    print deltaL[i], "   ", exp_from_deltaLto16[i], "   ", exp_from_deltaLto16_stdvs[i]
'''
