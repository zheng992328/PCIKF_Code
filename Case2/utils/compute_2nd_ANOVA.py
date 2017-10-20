# -*- coding: utf-8 -*-

import numpy as np
import os
from ANOVA_parallel import ANOVA_parallel
from write_values import write_ki
from read_values import read_obs,read_obs_alltime
from runexe import runexe
from para_key_mod import para_keys_modify
from MapParameters import MapParameters
from time_modify import time_modify
def compute_2nd_ANOVA(pair_all,iPair,C_O0_t,C_O1_t,A2,Pairs,iO,C_O0,C_O1,kl_term,out_node,Npair,C_O2,C_P0,C_P1,C_P2,ki_mean,root_directory,Nod_num,filename_KI,obs_Num,filename_result,filename_time,tt):
    newPair=1
    for i in xrange(Npair):
        if np.max(np.abs(pair_all[iPair,2:4]-Pairs[:,i].T))==0:
            newPair=0
            break
    CP1=np.zeros((kl_term,4))
    if newPair:
        Npair=Npair+1
        P=np.array([-np.sqrt(3),0,np.sqrt(3)])
        CP=np.zeros(kl_term)
        CP[int(pair_all[iPair,2])]=P[0]
        CP[int(pair_all[iPair,3])]=P[0]
        CP1[:,0]=CP
        
        CP=np.zeros(kl_term)
        CP[int(pair_all[iPair,2])]=P[2]
        CP[int(pair_all[iPair,3])]=P[0]
        CP1[:,1]=CP
        
        CP=np.zeros(kl_term)
        CP[int(pair_all[iPair,2])]=P[0]
        CP[int(pair_all[iPair,3])]=P[2]
        CP1[:,2]=CP
        
        CP=np.zeros(kl_term)
        CP[int(pair_all[iPair,2])]=P[2]
        CP[int(pair_all[iPair,3])]=P[2]
        CP1[:,3]=CP
        
        temp=np.zeros((4,out_node))
        for i in range(4):
#            y=np.squeeze(np.log(ki_mean)+np.dot(Fi[:,:kl_term],CP1[:,i]))
#            y=np.exp(y)
            p=MapParameters(CP1[:,i],kl_term,C_P0,C_P1,C_P2)
            np.savetxt(r'../para_updating_{0}_{1}_22.txt'.format(tt,i),p)
            p=np.exp(p)
            p=para_keys_modify(p)
            write_ki(Nod_num,p,i,filename_KI,root_directory)
            time_modify(tt,i,filename_time,root_directory)
            # if np.mod(i,4)==0:
            print '2_nd_{0}_{1}'.format(iO,i)
            runexe(i,root_directory)
            
            obs_p2=read_obs_alltime(tt,0,obs_Num,Nod_num,root_directory,filename_result)
            temp[i,:]=obs_p2-ANOVA_parallel(CP1[:,i],C_O0,C_O1,kl_term)
        
        b=np.zeros((9,out_node))
        b[0,:]=temp[0,:]
        b[2,:]=temp[1,:]
        b[6,:]=temp[2,:]
        b[8,:]=temp[3,:]
        
        k=np.dot(np.linalg.inv(A2),b)
        for s in range(out_node):
            t=k[:,s].reshape((3,3))
            C_O0_t[s]=C_O0_t[s]+t[0,0]
            C_O1_t[s,:,int(pair_all[iPair,3])]=C_O1_t[s,:,int(pair_all[iPair,3])]+t[0,1:3]
            C_O1_t[s,:,int(pair_all[iPair,2])]=C_O1_t[s,:,int(pair_all[iPair,2])]+t[1:3,0].T
            C_O2[s,:,:,int(pair_all[iPair,2]),int(pair_all[iPair,3])]=t[1:3,1:3]
            if s==iO:
                pair_all[iPair,1]=np.sum(t[1:3,1:3]*t[1:3,1:3])
        
        Ratio=np.sum(pair_all[:iPair,1])/(np.sum(pair_all[:iPair,0])+1e-10)
        Err_2nd=Ratio*np.sum(pair_all[iPair:,0])
        Pairs=np.hstack((Pairs,pair_all[iPair,2:4].reshape(-1,1)))  ##横向拼接，（2,1）与（2,1）的才可以拼接，与（2，）是拼接不起来的
    else:
        var=C_O2[iO,:,:,int(pair_all[iPair,2]),int(pair_all[iPair,3])]**2
        pair_all[iPair,1]=np.sum(var.reshape((4,1)))
        Ratio=np.sum(pair_all[:iPair,1])/(np.sum(pair_all[:iPair,0])+(1e-10))
        Err_2nd=Ratio*np.sum(pair_all[iPair:,0])
    return Err_2nd,C_O0_t,C_O1_t,Pairs,pair_all,Npair,C_O2
        


            
            
            
            