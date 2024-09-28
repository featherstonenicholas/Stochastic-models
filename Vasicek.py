import numpy as np
from scipy.stats import norm
#Price a ZCP using Vasicek model
# set up model parameters

theta=0.08
kappa=0.86
sigma=0.01
r_0=0.06
T_0=1
T_1=2

#price a bond using vasicek model

def vasicek_bond(t,T,theta,kappa,sigma,r_0):
    #calculate sub parts
    A=r_0*(1-np.exp(-kappa*(T-t)))/kappa
    B=(theta-sigma**2/(2*kappa**2))*( (1-np.exp(-kappa*(T-t))) /kappa -(T-t)  )
    C=(sigma**2/(4*kappa))*( 1-np.exp(-kappa*(T-t)))**2
    
    return np.exp(-A +B - C)

P_0=vasicek_bond(0,T_0,theta,kappa,sigma,r_0)
P_1=vasicek_bond(0,T_1,theta,kappa,sigma,r_0)


# calculare convexity adjustment T_0 to T_1

def gamma(t,T_0,T_1,theta,kappa,sigma,r_0):
    P_0=vasicek_bond(0,1,theta,kappa,sigma,r_0)
    P_1=vasicek_bond(0,5/4,theta,kappa,sigma,r_0)
    C=(sigma**2/kappa**3)*( 1-1.5*np.exp(-kappa*(T_1-T_0)) +0.5*np.exp(-2*kappa*(T_1-T_0))- np.exp(-kappa*(T_0-t))+np.exp(-kappa*(T_1-t))-0.5*np.exp(-2*kappa*(T_1-t)) +0.5*np.exp(-kappa*(T_1+T_0-2*t)))
    result=(P_0/P_1)*(np.exp(C)-1)*(1/(T_1-T_0))
    return result


#price a european call on ZCB using Vasicek. Price of t=0 call, maturity T_0 on a T_1 bond

#K=P_1/P_0

def vasicek_call(T_0,T_1,K,theta,kappa,sigma,r_0):
    P_0=vasicek_bond(0,T_0,theta,kappa,sigma,r_0)
    P_1=vasicek_bond(0,T_1,theta,kappa,sigma,r_0)
    
    vol=(sigma**2/kappa**2)*((np.exp(-kappa*T_0)-np.exp(-kappa*T_1))**2)*(np.exp(2*kappa*T_0)-1)/(2*kappa)
    d_1=(np.log(P_1/(K*P_0))+0.5*vol)/np.sqrt(vol)
    d_2=(np.log(P_1/(K*P_0))-0.5*vol)/np.sqrt(vol)
    
    call=P_1*norm.cdf(d_1) - K*P_0*norm.cdf(d_2)
   
    return call

ZCB=np.ones(121)
T0=0.25*np.arange(121)
delta=0.25
for i in range(121):
    ZCB[i]=vasicek_bond(0,T0[i],theta,kappa,sigma,r_0)

K = (ZCB[1] - ZCB[120]) / (delta*sum(ZCB[2:]))

def myIntegral(b,nu,t0,t1):
    return nu**2/b**2 *(np.exp(-b*t0) - np.exp(-b*t1))**2 * (np.exp(2*b*t0)-1)/(2*b)
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
#print(CapVasicek(kappa,sigma,K,T0,30,delta,ZCB))
vol=(sigma**2/kappa**2)*((np.exp(-kappa*T_0)-np.exp(-kappa*T_1))**2)*(np.exp(2*kappa*T_0)-1)/(2*kappa)
mu=np.log(P_0/P_1)+vol
print(np.sqrt(vol))