# -*- coding: utf-8 -*-
from HermiteP import HermiteP
import numpy as np
def ANOVA_parallel(CP,C0,C1,NR):
    '''
    CP: collocation points
    C0: 0th order coefficient
    C1: 1st order coefficient
    '''
    F1=C0
    for i in range(NR):
        t=HermiteP(CP[i])
        F1=F1+np.dot(C1[:,:,i],t[1:3])
    return F1