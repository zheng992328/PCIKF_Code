# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import os,multiprocessing
import matplotlib.pyplot as plt
from write_values import write_ki
from read_values import read_obs
from runexe import runexe
from HermiteP import HermiteP
from para_key_mod import para_keys_modify
from MapParameters import MapParameters
# from compute_2nd_ANOVA_1 import compute_2nd_ANOVA
from KLexpansion import KLexpansion
from time_modify import time_modify
from cal_ANOVA_11 import cal_ANOVA_1

def get_PCE_coef(kl_term,Nod_num,C_P0,C_P1,obs_Num,sigma,ki_mean,deltax,deltay,dx,dy,x,y,para_for_ogs,filename_result,filename_KI,filename_time,root_directory,out_node,tt):
    ##计算ANOVA中的常数项
    CP0=np.zeros((kl_term,1))
    p=MapParameters(CP0,kl_term,C_P0,C_P1)
    # np.savetxt(r'../para_updating_{0}_0.txt'.format(tt),p)
    p=np.exp(p)
    p=para_keys_modify(p)
    write_ki(Nod_num,p,0,filename_KI,root_directory)
    time_modify(tt,0,filename_time,root_directory)
    print 'start 0 order'
    runexe(0,root_directory)
    obs_p2=read_obs(Nod_num,obs_Num,0,root_directory,filename_result)
    C_O0=obs_p2
    
    
    ##计算ANOVA中的一阶项
    P=np.array([-np.sqrt(3),0,np.sqrt(3)])
    PC1=np.zeros((3,3))
    for i in range(3):
        PC1[:,i]=HermiteP(P[i])
    
    C_O1=np.zeros((out_node,2,kl_term))
    C_O0_a=np.zeros_like(C_O0)  #将一阶分解中的常数项并到之前求的常数项中去
    
    A1=PC1.T
    Energy0=np.zeros_like(C_O1)
    b=np.zeros((3*kl_term,out_node))
    # for i in xrange(kl_term):
        # b=np.zeros((3,out_node))
        # CP=np.zeros((kl_term,1))
        # CP[i,0]=P[0]
# #        y=np.squeeze(np.log(ki_mean)+np.dot(Fi[:,:kl_term],CP))
# #        y=np.exp(y)
        # p=MapParameters(CP,kl_term,C_P0,C_P1)
        # # np.savetxt(r'./para_updating_{0}_{1}_11.txt'.format(tt,i),p)
        # p=np.exp(p)
        # if np.mod(i,5)==0:
            # print 'i={0}_1'.format(i)
        # write_ki(Nod_num,p,0,filename_KI,root_directory)
        # time_modify(tt,0,filename_time,root_directory)
        # runexe(0,root_directory)
        
        # obs_p2=read_obs(Nod_num,obs_Num,0,root_directory,filename_result)
        # b[0,:]=obs_p2-C_O0
        
        # CP=np.zeros((kl_term,1))
        # CP[i,0]=P[2]
# #        y=np.squeeze(np.log(ki_mean)+np.dot(Fi[:,:kl_term],CP))
# #        y=np.exp(y)
        # p=MapParameters(CP,kl_term,C_P0,C_P1)
        # # np.savetxt(r'./para_updating_{0}_{1}_12.txt'.format(tt,i),p)
        # p=np.exp(p)
        # if np.mod(i,5)==0:
            # print 'i={0}_2'.format(i)
        # write_ki(Nod_num,p,0,filename_KI,root_directory)
        # time_modify(tt,0,filename_time,root_directory)
        # runexe(0,root_directory)
        
        # obs_p2=read_obs(Nod_num,obs_Num,0,root_directory,filename_result)
        # b[2,:]=obs_p2-C_O0
        
        # k=np.dot(np.linalg.inv(A1),b)
        # C_O1[:,:,i]=k[1:3,:].T
        # C_O0_a=C_O0_a+k[0,:].T
        # Energy0[:,:,i]=C_O1[:,:,i]*C_O1[:,:,i]
        # # bb=np.vstack((bb,b))
    # np.savetxt('bb.txt',bb)
    print 'start 1 order'
    pool=multiprocessing.Pool(20)
    result=[]
    for i in xrange(kl_term):
        result.append(pool.apply_async(cal_ANOVA_1,(i,kl_term,out_node,Nod_num,obs_Num,C_P0,C_P1,C_O0,filename_KI,filename_time,filename_result,root_directory,tt)))
    pool.close()
    pool.join()
    for i in xrange(len(result)):
        b[3*i:3*i+3,:]=result[i].get()
    # np.savetxt('b.txt',b)
    for i in xrange(kl_term):
        k=np.dot(np.linalg.inv(A1),b[3*i:3*i+3,:])
        C_O1[:,:,i]=k[1:3,:].T
        C_O0_a=C_O0_a+k[0,:].T
        Energy0[:,:,i]=C_O1[:,:,i]+C_O1[:,:,i]
    C_O0=C_O0+C_O0_a
    print 'over'
    
    # ##计算ANOVA中的二阶项
    # C_O2=np.zeros((out_node,2,2,kl_term,kl_term))
    # C_O0_t=C_O0
    # C_O1_t=C_O1
    
    # A2=np.zeros((9,9))
    # n1=0
    # for u in range(3):
        # for v in range(3):
            # t=np.dot(PC1[:,v].reshape((-1,1)),PC1[:,u].reshape((1,-1)))
            # A2[n1,:]=t.reshape((1,9))
            # n1=n1+1
    
    
    # #二阶ANOVA自适应标准
    # Pairs=np.zeros((2,0))
    # Npair=0
# #    Energy_1st_t=np.sum(np.sum(Energy0,axis=1),axis=1)
    # Energy_1st=np.sum(Energy0,axis=1)
    # N_Total_pair=kl_term*(kl_term-1)/2
    # Thresh=0.5
    # for iO in xrange(out_node):
        # pair_all=np.zeros((N_Total_pair,4))
        # n1=0
        # for u in range(kl_term-1):
            # for v in range(u+1,kl_term):
                # pair_all[n1,:]=np.array([Energy_1st[iO,u]*Energy_1st[iO,v],0,u,v])
        # pair_all=pd.DataFrame(pair_all,columns=['a','b','c','d'],index=[ind for ind in range(N_Total_pair)])
        # pair_all=pair_all.sort_index(by='a',ascending=False)  ##按第一列降序排列
        # pair_all=pair_all.values
        
        # MaxPair=N_Total_pair
        # for iPair in range(MaxPair):
            # Err_2nd,C_O0_t,C_O1_t,Pairs,pair_all,Npair,C_O2=compute_2nd_ANOVA(pair_all,iPair,C_O0_t,C_O1_t,A2,Pairs,iO,C_O0,C_O1,kl_term,out_node,Npair,C_O2,C_P0,C_P1,C_P2,ki_mean,root_directory,Nod_num,filename_KI,obs_Num,filename_result,filename_time,tt)
            # Total_out_E=np.sum(C_O1_t[iO,:,:]**2)+np.sum(C_O2[iO,:,:,:,:]**2)  #sum如果不指定axis则默认所有元素相加
# #            str_='{0} {1} {2}'.format(iO,MaxPair,iPair)
# #            print str_
            # if Err_2nd<Thresh*(Total_out_E+0.05):
                # break
    # print '2nd ANOVA is over'
    # C_O0=C_O0_t
    # C_O1=C_O1_t
    return C_O0,C_O1

if __name__=='__main__':
    Obs_Num=16
    Nod_num=1681
    obs_Num=[213,221,229,237,623,631,639,647,1033,1041,1049,1057,1443,1451,1459,1467]
    time_step=20
    out_node=Obs_Num
    kl_term=25
    sigma=0.6
    ki_mean=2.9e-13
    x=2.0*30
    y=1.0*30
    deltax=x/3.0
    deltay=y/3.0
    dx=0.05*30
    dy=0.025*30
    test_num=20
    m=int(y/dy+1)
    n=int(x/dx+1)
    para_for_ogs='water_gas'
    filename_result='water_gas_domain_quad.tec'
    filename_KI='water_gas_Ki.direct'
    filename_time='water_gas.tim'
    root_directory=os.getcwd()
    root_directory=os.path.dirname(root_directory)
    
    
    # def cum(Fi,m):
        # sum_numta=0
        # for i in xrange(m):
            # sum_numta+=Fi[i,i]
        # return sum_numta
      
    Fi=KLexpansion(m,n,x,y,deltax,deltay,sigma)
    # reserved=x*y*(sigma**2)*0.8
    # for m in xrange(Fi.shape[1]):
        # if cum(Fi,m)-reserved>0:
            # print m
            # break
    
    
    
    
    Y_mean=np.ones((Nod_num,1))*np.log(ki_mean)
    C_P0=Y_mean
    C_P1=np.zeros((Nod_num,2,kl_term))
    C_P1[:,0,:]=Fi[:,:kl_term]
    # C_P2=np.zeros((Nod_num,2,2,kl_term,kl_term))
    corrcoef=[]
    for t in xrange(1,time_step+1):
        C_O0,C_O1=get_PCE_coef(kl_term,Nod_num,C_P0,C_P1,obs_Num,sigma,ki_mean,deltax,deltay,dx,dy,x,y,para_for_ogs,filename_result,filename_KI,filename_time,root_directory,out_node,t)
        
        yy1=np.zeros((out_node,test_num))
        yy2=np.zeros((out_node,test_num))
        k1=np.random.normal(0,1,(kl_term,test_num))
        for k in xrange(test_num):
            # print 'k={0}'.format(k)
            k1_1=np.squeeze(k1[:,k])
            yy1[:,k]=MapParameters(k1_1,kl_term,C_O0,C_O1)
          
        y=np.log(ki_mean)+np.dot(Fi[:,:kl_term],k1)
        for i in xrange(test_num):
            # print i
            p=np.exp(y[:,i])
            p=para_keys_modify(p)
            write_ki(Nod_num,p,i,filename_KI,root_directory)
        jobs=[]
        for i in xrange(test_num):
            pp=multiprocessing.Process(target=runexe,args=(i,root_directory))
            jobs.append(pp)
            pp.start()
        for j in jobs:
            j.join()
        
        for i in xrange(test_num):
            yy2[:,i]=read_obs(Nod_num,obs_Num,i,root_directory,filename_result)
            
        result_1=yy1.reshape((out_node*test_num,1))
        result_2=yy2.reshape((out_node*test_num,1))
        corrcoef_=np.corrcoef(np.squeeze(result_1),np.squeeze(result_2))
        corrcoef.append(corrcoef_[0,1])
    plt.plot(range(1,time_step+1),corrcoef)
        # print np.corrcoef(np.squeeze(result_1),np.squeeze(result_2))
    # plt.scatter(result_1,result_2,marker='o')
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    