# -*- coding: utf-8 -*-

import numpy as np
from runexe import runexe
from HermiteP import HermiteP
from MapParameters import MapParameters
from write_values import write_ki
from read_values import read_obs
from para_key_mod import para_keys_modify
from time_modify import time_modify

def cal_ANOVA_1(i,kl_term,out_node,Nod_num,obs_Num,C_P0,C_P1,C_O0,filename_KI,filename_time,filename_result,root_directory,tt):
    P=np.array([-np.sqrt(3),0,np.sqrt(3)])
    b=np.zeros((3,out_node))
    CP=np.zeros((kl_term,1))
    CP[i,0]=P[0]
    p=MapParameters(CP,kl_term,C_P0,C_P1)
    # np.savetxt('para_updating_{0}_11.txt'.format(tt),p)
    p=np.exp(p)
    p=para_keys_modify(p)
    write_ki(Nod_num,p,i,filename_KI,root_directory)
    time_modify(tt,i,filename_time,root_directory)
    runexe(i,root_directory)
    obs_p2=read_obs(Nod_num,obs_Num,i,root_directory,filename_result)
    b[0,:]=obs_p2-C_O0
    
    CP=np.zeros((kl_term,1))
    CP[i,0]=P[2]
    p=MapParameters(CP,kl_term,C_P0,C_P1)
    # np.savetxt('para_updating_{0}_12.txt'.format(tt),p)
    p=np.exp(p)
    p=para_keys_modify(p)
    write_ki(Nod_num,p,i,filename_KI,root_directory)
    time_modify(tt,i,filename_time,root_directory)
    runexe(i,root_directory)
    
    obs_p2=read_obs(Nod_num,obs_Num,i,root_directory,filename_result)
    b[2,:]=obs_p2-C_O0
    return b
    
    


