# -*- coding: utf-8 -*-
import numpy as np
from Findwxz import Findwxz
import time
import matplotlib.pyplot as plt
def KLexpansion(m,n,Lx,Ly,CorLengthx,CorLengthy,sigma):
    """
    Lx,Ly are length in x,y direction respectively\n
    m,n are rows and columns\n
    CorLenthx,CorLengthy are correlation length\n
    sigma is stadard deviation
    """


    Nrow=m
    Ncol=n
    Npd=500  ##the number in fun.dat
    Nroot=40
    dx=Lx/(Ncol-1)*np.ones(Ncol)
    dy=Ly/(Nrow-1)*np.ones(Nrow)
    
    X0RLn=0.1
    dXRLn=1.0e-5
    z0RLn=0.1
    dzRLn=1.0e-5
    Rx=CorLengthx/Lx
    Ry=CorLengthy/Ly

    wx=Findwxz(Nroot,X0RLn,dXRLn,1.0,Rx)
    wy=Findwxz(Nroot,z0RLn,dzRLn,1.0,Ry)

    for i in xrange(Nroot):
        wx[i]=wx[i]/Lx
        wy[i]=wy[i]/Ly

    RLn=np.zeros((3,Nroot**2))
    for i in xrange(Nroot):
        for j in xrange(Nroot):
            k=i*Nroot+j   
            RLn[0,k]=wx[i]
            RLn[1,k]=wy[j]
            RLn[2,k]=4.0*CorLengthx*CorLengthy*sigma/((CorLengthx*wx[i])**2+1.0)/((CorLengthy*wy[j])**2+1.0)
            if k>0:
                for kk in xrange(k,0,-1):
                    if RLn[2,kk]>RLn[2,kk-1]:
                        r1=RLn[0,kk]
                        r2=RLn[1,kk]
                        r3=RLn[2,kk]
                        RLn[0,kk]=RLn[0,kk-1]
                        RLn[1,kk]=RLn[1,kk-1]
                        RLn[2,kk]=RLn[2,kk-1]
                        RLn[0,kk-1]=r1
                        RLn[1,kk-1]=r2
                        RLn[2,kk-1]=r3
                        

    Fi=np.zeros((Npd,Ncol*Nrow))
    for k in xrange(Npd):
        y=0
        for j in xrange(Nrow):
            x=0
            for i in xrange(Ncol):
                N=j*Ncol+i
                r1=CorLengthx*RLn[0,k]
                r2=CorLengthy*RLn[1,k]
                r4=RLn[0,k]*x
                r5=RLn[1,k]*y
                r3=(r1*np.cos(r4)+np.sin(r4))*(r2*np.cos(r5)+np.sin(r5))/np.sqrt(((r1*r1+1.0)*Lx/2.0+CorLengthx)*((r2*r2+1.0)*Ly/2.0+CorLengthy))
                Fi[k,N]=np.sqrt(RLn[2,k])*r3
                x=x+dx[i]
            y=y+dy[j]
    Eigenv=np.zeros(Npd)
    for i in xrange(Npd):
        Eigenv[i]=RLn[2,i]
    Fi=Fi.T
	###要测试取多少kl项合适的时候把下面注释取消
    # cum=np.zeros_like(Eigenv)
    # c=np.zeros_like(cum)
    # for i in xrange(1,Npd+1):
        # cum[i-1]=np.sum(Eigenv[:i])
    # c=cum/cum[-1]
    # plt.plot(range(Npd),c)
    # plt.show()
    return Fi

if __name__=='__main__':
    m=n=41
    Lx=30.0
    Ly=15.0
    CorLengthx=Lx/3
    CorLengthy=Ly/3
    sigma=0.6
    Fi=KLexpansion(m,n,Lx,Ly,CorLengthx,CorLengthy,sigma)
                
    
    
    
    