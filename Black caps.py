import numpy as np
from scipy.stats import *
from scipy.optimize import *

ForwardRates = [0.06,0.08,0.09,0.10,0.10,0.10,0.09,0.09] 
delta=0.25
T0 = np.arange(0,2.25,delta) 

Maturities = [2]
CapPrices = [0.01]
#create the ZCB prices from forward rates

def discount(times,ForwardRates,time_delta):
    delta=time_delta
    bondprices=np.ones(len(times))
    for i in range(1, len(times)):
        bondprices[i] = bondprices[i-1]/(1 + delta*ForwardRates[i-1])
    return bondprices

ZeroBondPrices=discount(T0,ForwardRates,delta)

# Create an array withh the 6 month forward swap rates

def forward_swap(ZeroBondPrices,delta,Maturities):
    kappa=np.zeros(len(Maturities))
    
    for i in range(len(Maturities)):
        m=Maturities[i]
        jump=int(1/delta)
        kappa[i] = (ZeroBondPrices[1] - ZeroBondPrices[jump*m]) / (delta*sum(ZeroBondPrices[2:(jump*m+1)]))
        
    return kappa

kappa=forward_swap(ZeroBondPrices,delta,Maturities)

#function to calculate sigma for vasicek HJM
def myIntegral(b,nu,t0,t1):
    return nu**2/b**2 *(np.exp(-b*t0) - np.exp(-b*t1))**2 * (np.exp(2*b*t0)-1)/(2*b)

#M is maturity of cap, T is time vector, Z is discount vector, k is the par rate for cap
def CapVasicek(b,nu,k,T,M,delta,Z):
    myAns = 0
    jump=int(1/delta)
    for i in range(1, (jump*M)):
        I = myIntegral(b, nu, T[i], T[i+1])
        
        d1 = (np.log(Z[i+1]/Z[i]*(1+delta*k)) + 0.5*I)/np.sqrt(I)
        d2 = (np.log(Z[i+1]/Z[i]*(1+delta*k)) - 0.5*I)/np.sqrt(I)
        cplt_i = Z[i]*norm.cdf(-d2,0,1)-(1+delta*k)*Z[i+1]*(norm.cdf(-d1,0,1))
        #print(Z[i+2])
        #print(cplt_i)
        myAns = myAns + cplt_i 
    return myAns

def BlackCap(sig,k,fwds,T,M,delt,Z):
    myAns = 0
    
    jump=int(1/delta)
    
    for i in range(1, (jump*M)):
        d1 = (np.log(fwds[i]/k) + 0.5*(sig**2)*T[i])/(sig*np.sqrt(T[i]))
        d2 = d1 - sig*np.sqrt(T[i])
        
        cplt_i = delt*Z[i+1]*(fwds[i]*norm.cdf(d1,0,1) - k*norm.cdf(d2,0,1))
        #print('i, fwds[i], k, T[i-1], T[i], sig', i, fwds[i], k, T[i-1], T[i], sig)
        #print('d1, d2, cplt_i', d1, d2, cplt_i)
        myAns = myAns + cplt_i
    return myAns

def BlackVega(delt,Z,fwds,T,sig,M,k):
    myAns = 0
    jump=int(1/delta)
    for i in range(1, (jump*M)):
        #print(fwds[i])
        #print(sig)
        d1 = (np.log(fwds[i]/k) + 0.5*sig**2*T[i])/(sig*np.sqrt(T[i]))
        #d2 = d1 - sig*sqrt(T[i])
        cplt_Vega = delt*Z[i+1]*fwds[i]*np.sqrt(T[i])*norm.pdf(d1,0,1)
        myAns = myAns + cplt_Vega
    return myAns

def black_imp_vol(ForwardRates,T0,m,delta,kappa,ZeroBondPrices,CapPrice):
    BCap = lambda iv: BlackCap(iv,kappa,ForwardRates,T0,m,delta,ZeroBondPrices)-CapPrice
    imp=bisect(BCap,0.005,10)
    return imp
print('Cap Prices',CapPrices)
print('Maturities', Maturities)
print('Forward rates',ForwardRates)
print('Discounts',ZeroBondPrices)
print('time',T0)
print('Forward swaps',kappa)


print(BlackCap(0.141,kappa[0],ForwardRates,T0,2,delta,ZeroBondPrices))
print(black_imp_vol(ForwardRates,T0,2,delta,kappa[0],ZeroBondPrices,CapPrices[0]))