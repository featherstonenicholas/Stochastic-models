import numpy as np
from scipy.stats import norm
#Price a ZCP using Vasicek model
# set up model parameters

theta=0.09
kappa=0.86
sigma=0.0148
r_0=0.08
T_0=1
T_1=1.25

#price a bond using vasicek model

def vasicek_bond(t,T,theta,kappa,sigma,r_0):
    #calculate sub parts
    B=(1-np.exp(-kappa*(T-t)))/kappa
    A_1=(theta-sigma**2/(2*kappa**2))*( B - (T-t))
    A_2=(sigma**2/(4*kappa))*B**2
    
    return np.exp(-r_0*B+A_1+A_2)

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

K=P_1/P_0

def vasicek_call(T_0,T_1,K,theta,kappa,sigma,r_0):
    P_0=vasicek_bond(0,T_0,theta,kappa,sigma,r_0)
    P_1=vasicek_bond(0,T_1,theta,kappa,sigma,r_0)
    
    vol=(sigma**2/kappa**2)*((np.exp(-kappa*T_0)-np.exp(-kappa*T_1))**2)*(np.exp(2*kappa*T_0)-1)/(2*kappa)
    d_1=(np.log(P_1/(K*P_0))+0.5*vol)/np.sqrt(vol)
    d_2=(np.log(P_1/(K*P_0))-0.5*vol)/np.sqrt(vol)
    
    call=P_1*norm.cdf(d_1) - K*P_0*norm.cdf(d_2)
   
    return call
ZCB=np.ones(120)
T0=0.25*np.arange(120)
delta=0.25
for i in range(120):
    ZCB[i]=vasicek_bond(0,T0[i],theta,kappa,sigma,r_0)

K = (ZCB[1] - ZCB[2*m]) / (delta*sum(ZCB[2:]))
print(K)