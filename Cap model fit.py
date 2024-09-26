import numpy as np
from scipy.stats import *
from scipy.optimize import minimize, bisect

def myIntegral(b,nu,t0,t1):
    return nu**2/b**2 *(np.exp(-b*t0) - np.exp(-b*t1))**2 * (np.exp(2*b*t0)-1)/(2*b)

def CapVasicek(b,nu,kappa,ForwardRates,T,M,delta,Z):
    myAns = 0
    k = kappa[M]
    for i in range(1, (2*M+2)):
        I = myIntegral(b, nu, T[i], T[i+1])
        print(I)
        d1 = (np.log(Z[i+1]/Z[i]*(1+delta*k)) + 0.5*I)/np.sqrt(I)
        d2 = (np.log(Z[i+1]/Z[i]*(1+delta*k)) - 0.5*I)/np.sqrt(I)
        cplt_i = Z[i]*norm.cdf(-d2,0,1)-(1+delta*k)*Z[i+1]*(norm.cdf(-d1,0,1))
        #print(Z[i+2])
        #print(cplt_i)
        myAns = myAns + cplt_i 
    return myAns

def BlackVega(delt,Z,fwds,T,sig,M,kap):
    myAns = 0
    k = kap[M]
    for i in range(1, (2*M+2)):
        #print(fwds[i])
        #print(sig)
        d1 = (np.log(fwds[i]/k) + 0.5*sig**2*T[i])/(sig*np.sqrt(T[i]))
        #d2 = d1 - sig*sqrt(T[i])
        cplt_Vega = delt*Z[i+1]*fwds[i]*np.sqrt(T[i])*norm.pdf(d1,0,1)
        myAns = myAns + cplt_Vega
    return myAns

def BlackCap(sig,kap,fwds,T,M,delt,Z):
    myAns = 0
    k = kap[M]
    for i in range(1, (2*M+2)):
        d1 = (np.log(fwds[i]/k) + 0.5*(sig**2)*T[i])/(sig*np.sqrt(T[i]))
        d2 = d1 - sig*np.sqrt(T[i])
        cplt_i = delt*Z[i+1]*(fwds[i]*norm.cdf(d1,0,1) - k*norm.cdf(d2,0,1))
        #print('i, fwds[i], k, T[i-1], T[i], sig', i, fwds[i], k, T[i-1], T[i], sig)
        #print('d1, d2, cplt_i', d1, d2, cplt_i)
        myAns = myAns + cplt_i
    return myAns
    

# In[2]:

ForwardRates = [0.06,0.08,0.09,0.10,0.10,0.10,0.09,0.09] # c(0.06,0.08,0.09,0.10,0.10,0.10,0.09,0.09)
T0 = np.arange(0,4.5,0.5) #<- seq(0,4,by=0.5)

CapPrices = [0.002,0.008,0.012,0.016]
CapMat = [1,2,3,4]
delta = 0.5

ZeroBondPrices = [1.0]*len(T0)
d1 = [1.0]*8

for i in range(1, len(ZeroBondPrices)):
    ZeroBondPrices[i] = ZeroBondPrices[i-1]/(1 + delta*ForwardRates[i-1])

Maturities = [1,2,3,4]
kappa = [0,0,0,0]
for m in Maturities:
    #print(m) # m=1, T0=1/2, Tn=m
    kappa[m-1] = (ZeroBondPrices[1] - ZeroBondPrices[2*m]) / (delta*sum(ZeroBondPrices[2:(2*m+1)]))

#print(ForwardRates,T0,m,delta,ZeroBondPrices,CapPrices, kappa)

# In[3]:

impVol = []

for m in range(0,4):
    BCap = lambda iv: BlackCap(iv,kappa,ForwardRates,T0,m,delta,ZeroBondPrices)-CapPrices[m]
    impVol.append(bisect(BCap,0.005,0.25))
    #print(impVol)

print(impVol)

# In[4]:

Vegas = []
for m in range(0,4):
    Vegas.append(BlackVega(delta,ZeroBondPrices,ForwardRates,T0,impVol[m],m,kappa))
    #print(Vegas)
#print(Vegas)

# In[5]:

def objectiveFunction(par):
    b = par[0]
    nu = par[1]
    myAns = 0
    for m in range(4):
        Cn_model = CapVasicek(b,nu,kappa,ForwardRates,T0,m,delta,ZeroBondPrices)
        myAns = myAns + (1/Vegas[m]**2)*(Cn_model-CapPrices[m])**2
    return myAns


# In[6]:

#print(kappa,ForwardRates,T0,ZeroBondPrices)
#print(CapVasicek(1.03,0.03,kappa,ForwardRates,T0,3,delta,ZeroBondPrices))


# In[7]:

x = minimize(objectiveFunction, x0=[0.01,0.03])
#print(x)
print(kappa)
print(len(T0),len(ForwardRates))
