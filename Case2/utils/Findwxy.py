# -*- coding: utf-8 -*-

import numpy as np
import time

def fout(x,C,L):
    f=(np.power(C,2)*np.power(x,2)-1)*np.sin(x*L)-2*C*x*np.cos(x*L)
    return f
    
def Findwxz(Nroot,x0,dx,L,CorLength):
    accuricywn=1e-5
    Nr=0
    x1=x0
    f0=fout(x0,CorLength,L)
    x=np.zeros(Nroot)
    y=np.zeros(Nroot)
    while Nr<Nroot:
        x1=x1+dx
        f1=fout(x1,CorLength,L)
        if f1==0:
            x[Nr]=x1
            y[Nr]=fout(x[Nr],CorLength,L)
            Nr+=1
            x0=x1+0.5*dx
            f0=fout(x0,CorLength,L)
        elif f0*f1<0:
            x3=x1
            while np.abs(x0-x1)>accuricywn:
                x2=(x0+x1)*0.5
                f2=fout(x2,CorLength,L)
                if f2==0:
                    x[Nr]=x2
                    y[Nr]=f2
                    Nr+=1
                    x0=x1+0.5*dx
                    f0=fout(x0,CorLength,L)
                elif f1*f2<0:
                    x0=x2
                    f0=f2
                else:
                    x1=x2
                    f1=f2
            
            x[Nr]=x1-0.5*np.abs(x1-x0)
            y[Nr]=fout(x[Nr],CorLength,L)
            Nr+=1
            x0=x3
            f0=fout(x3,CorLength,L)
            
    return x


if __name__=='__main__':
    Nroot=40
    X0RLn=0.1
    dXRLn=1.0e-5
    z0RLn=0.1
    dzRLn=1.0e-5
    CorLengthx=100.0
    CorLengthy=50.0
    Lx=Ly=200.0
    Rx=CorLengthx/Lx
    Ry=CorLengthy/Ly
    start=time.time()
    x=Findwxz(Nroot,X0RLn,dXRLn,1.0,Rx)
    y=Findwxz(Nroot,X0RLn,dXRLn,1.0,Ry)
    print time.time()-start
            