# -*- coding: utf-8 -*-

import numpy as np
from HermiteP import HermiteP
def MapParameters(si,NR,C_P0,C_P1):
    '''
    mapping from random space to parameter space,using PCE representation of the parameters\n
    返回结果与C_P0有同样的维度，如果C_P0是测试PCE时候用的模型预测值，那返回的就是PCE求得的模型预测值
    如果C_P0是参数的先验，那返回的仍是参数
    '''
    # f=lambda x:np.matrix(x)
    p=np.squeeze(C_P0)
    HP=np.zeros((3,NR))
    for i in range(NR):
        HP[:,i]=HermiteP(si[i])
        p=p+np.dot(C_P1[:,:,i],HP[1:3,i])
    # p,HP,C_P2=map(f,[p,HP,C_P2])
    # for i in range(NR):
        # for j in range(NR):
            # for u in range(2):
                # for v in range(2):
                    # temp=np.dot(C_P2[:,v,u,j,i],HP[1+v,j])
                    # p=p+np.dot(temp,HP[1+u,i])
    p=np.array(p)
    return p

if __name__=='__main__':
    from KLexpansion import KLexpansion
    Nod_num=1681
    ki_mean=2.9e-15
    kl_term=20
    kl_term=20
    sigma=0.6
    ki_mean=2.9e-15
    deltax=1.0
    deltay=0.5
    dx=0.05
    dy=0.025
    x=2.0
    y=1.0
    m=int(y/dy+1)
    n=int(x/dx+1)
    Y_mean=np.ones((Nod_num,1))*np.log(ki_mean)
    Fi=KLexpansion(m,n,x,y,deltax,deltay,sigma)
    C_P0=Y_mean
    C_P1=np.zeros((Nod_num,2,kl_term))
    C_P1[:,0,:]=Fi[:,:kl_term]
    C_P2=np.zeros((Nod_num,2,2,kl_term,kl_term))
    CP0=np.zeros((kl_term,1))
    p=MapParameters(CP0,kl_term,C_P0,C_P1,C_P2)