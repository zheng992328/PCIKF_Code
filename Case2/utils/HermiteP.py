# -*- coding: utf-8 -*-

import numpy as np
def HermiteP(x):
    y=np.zeros(3)
    y[0]=1
    y[1]=x
    y[2]=(x**2-1)/np.sqrt(2)
    return y