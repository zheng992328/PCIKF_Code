# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
import time
import subprocess
import multiprocessing
import os,sys
import matplotlib.pyplot as plt
from utils.generate_obs_ref import generate_obs,write_KI
from utils.para_key_mod import para_keys_modify
from utils.write_values import write_ki
from utils.read_values import read_obs_alltime
from utils.read_values import read_obs
from utils.time_modify import time_modify
from utils.KLexpansion import KLexpansion
from utils.runexe import runexe

N=np.int(sys.argv[1])
index=np.int(sys.argv[2])
start=time.ctime()
time_step=20
Nod_num=1681
varR=80.0
obs_Num=[213,221,229,237,623,631,639,647,1033,1041,1049,1057,1443,1451,1459,1467]
obs_num=len(obs_Num)
sigma=0.6
ki_mean=2.9e-13
x=2.0*30
y=1.0*30
deltax=x/3.0
deltay=y/3.0
dx=0.05*30
dy=0.025*30
m=int(y/dy+1)
n=int(x/dx+1)
kl_term=25
root_directory=os.getcwd()
root_directory_true_obs=os.path.join(root_directory,'true_obs')
para_for_ogs='water_gas'
filename_result='water_gas_domain_quad.tec'
filename_KI='water_gas_Ki.direct'
filename_time='water_gas.tim'

#filename_Pressure='gas_PRESSURE1.direct'
beta=0.05
eps1=0.03
eps2=1e-2

## 产生观测值
Y_mean=np.ones((Nod_num,1))*np.log(ki_mean)
kl_true=np.random.standard_normal((kl_term,1))
Fi=KLexpansion(m,n,x,y,deltax,deltay,sigma)
Y_true=Y_mean+np.dot(Fi[:,:kl_term],kl_true)
np.savetxt('para_true.txt',Y_true)

ki_true=np.exp(Y_true)
ki_true=para_keys_modify(ki_true)
write_ki(Nod_num,ki_true,0,filename_KI,root_directory)
time_modify(time_step,0,filename_time,root_directory)
runexe(0,root_directory)
print time.ctime()
print 'get obs over'
obs=read_obs_alltime(time_step,0,obs_Num,Nod_num,root_directory,filename_result)

##调取观测值样本
obs_Pressure2=generate_obs(time_step,obs_num,varR,N,obs)
t=1
y_p1=obs_Pressure2[t]   ##尝试只用气压观测数据，没有基质吸力监测数据。
y_obs_p1=pd.DataFrame(y_p1)
y_obs_p1=y_obs_p1.T
y_obs_p1=y_obs_p1.values
y_obs_p1_ave=np.mean(y_obs_p1,axis=1)

#产生初始参数
kl_initial=np.random.standard_normal((kl_term,N))   ###是将它估计准了，就可以间接估计准整个参数场吗？？？
Fi=KLexpansion(m,n,x,y,deltax,deltay,sigma)
para_initial=np.log(ki_mean)+np.dot(Fi[:,:kl_term],kl_initial)
# para_initial=para_init(sigma,deltax,deltay,dx,dy,m,n,Nod_num,N,ki_mean)
kl_initial_ave=np.mean(kl_initial,axis=1)


#产生初始p1
p1_initial=np.zeros([Nod_num,N])
p1_initial=p1_initial+101325

##获取误差矩阵
y_obs_p1_error=np.zeros_like(y_obs_p1)
kl_error=np.zeros_like(kl_initial)
for i in xrange(N):
    y_obs_p1_error[:,i]=y_obs_p1[:,i]-y_obs_p1_ave[:]
    kl_error[:,i]=kl_initial[:,i]-kl_initial_ave[:]

Cd=np.dot(y_obs_p1_error,y_obs_p1_error.T)/(N-1)
Cm=np.dot(kl_error,kl_error.T)/(N-1)
mpr=kl_initial
m1=kl_initial
np.savetxt('para_initial.txt',kl_initial)
Imax=1

y_obs_prediction=np.zeros_like(y_obs_p1)
y_obs_prediction_error=np.zeros_like(y_obs_prediction)
p_after=np.zeros([Nod_num,N])
parY=np.zeros((kl_term,time_step))
updated_para=np.zeros((Nod_num,time_step))
spread=np.zeros(time_step)
rmse=np.zeros(time_step)
iter_num=[]
time_sum=np.sum(np.arange(1,21))

    
while True:
    print '1_time_{0}_iteration'.format(Imax)
    mat_to_arr=lambda x:np.array(x)
    m1,y_obs_prediction,y_obs_p1=map(mat_to_arr,[m1,y_obs_prediction,y_obs_p1])
    m1_average=np.mean(m1,axis=1)
    m1_error=np.zeros_like(m1)
    for i in xrange(N):
        m1_error[:,i]=m1[:,i]-m1_average[:]
        
    ks=np.exp(np.log(ki_mean)+np.dot(Fi[:,:kl_term],m1))
    for i in xrange(N):
        time_modify(1,i,filename_time,root_directory)
        ki_modify=ks[:,i]
        ki_modify=para_keys_modify(ki_modify)
        write_ki(Nod_num,ki_modify,i,filename_KI,root_directory) 
        

    jobs=[]
    for i in xrange(N):
        p=multiprocessing.Process(target=runexe,args=(i,root_directory,))
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()
   
    for k in xrange(N):
        Obs_p1=read_obs(Nod_num,obs_Num,k,root_directory,filename_result)
        for ii in xrange(obs_num):
            y_obs_prediction[ii][k]=Obs_p1[ii]  #组成观测点处的预测值矩阵
       
    # np.savetxt('y_obs_prediction_1.txt',y_obs_prediction)    
    y_obs_prediction_average=np.mean(y_obs_prediction,axis=1)
    for i in xrange(N):
        y_obs_prediction_error[:,i]=y_obs_prediction[:,i]-y_obs_prediction_average[:]
    
    wn_M=np.linalg.pinv(m1_error)
    G=np.dot(y_obs_prediction_error,wn_M)
    arr_to_mar=lambda x: np.matrix(x)
    mpr,m1,Cm,G,Cd,y_obs_prediction,y_obs_p1=map(arr_to_mar,[mpr,m1,Cm,G,Cd,y_obs_prediction,y_obs_p1])
    m2=beta*mpr+(1-beta)*m1-beta*Cm*(G.T)*(np.linalg.inv(Cd+G*Cm*G.T))*(y_obs_prediction-y_obs_p1-G*(m1-mpr))
    S1=np.trace(((y_obs_prediction-y_obs_p1).T)*np.linalg.inv(Cd)*(y_obs_prediction-y_obs_p1))
    m2=np.array(m2)
    y_obs_prediction=np.array(y_obs_prediction)
    # np.savetxt('m2.txt',m2)

   
    ks2=np.exp(np.log(ki_mean)+np.dot(Fi[:,:kl_term],m2))
 
    for i in xrange(N):
        ki_modify=ks2[:,i]
        ki_modify=para_keys_modify(ki_modify)
        write_ki(Nod_num,ki_modify,i,filename_KI,root_directory)  #para为初始参数，样本之间存在随机扰动


    jobs=[]
    for i in xrange(N):
        p=multiprocessing.Process(target=runexe,args=(i,root_directory,))
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()
   
    for k in xrange(N):
        Obs_p1=read_obs(Nod_num,obs_Num,k,root_directory,filename_result)
        for ii in xrange(obs_num):
            y_obs_prediction[ii][k]=Obs_p1[ii]  #组成观测点处的预测值矩阵
    # np.savetxt('y_obs_prediction_2.txt',y_obs_prediction)
    y_obs_prediction=np.matrix(y_obs_prediction)
 
    S2=np.trace(((y_obs_prediction-y_obs_p1).T)*np.linalg.inv(Cd)*(y_obs_prediction-y_obs_p1))
    print 'S1:',S1
    print 'S2:',S2
    print 'max(m2-m1):',np.max(np.abs(m2-m1))
    print 'S2-S1:',S2-S1
    print 'eps2*S1',eps2*S1,'\n'
    
    if S2<S1:
        if (np.max(np.abs(m2-m1))<eps1 or S1-S2<eps2*S1 or Imax>3):
            break
        m1=m2
        beta=2*beta
        if beta>1:
            beta=1
    else:
        beta=0.5*beta
    Imax+=1
    if (Imax>5):
        break
iter_num.append(Imax*1.0/time_sum)
parY[:,0]=np.average(m2,axis=1)
updated_para_each_time=np.log(ki_mean)+np.dot(Fi[:,:kl_term],m2)
spread[0]=np.sqrt(np.average(np.var(updated_para_each_time,axis=1)))
updated_para_average=np.mean(updated_para_each_time,axis=1)
updated_para[:,0]=updated_para_average
rmse[0]=np.sqrt(np.sum((updated_para_average-np.squeeze(Y_true))**2)/(Nod_num-1))

###这里得到了当前时间步的更新参数，要用它去计算目标函数
for t in range(2,time_step+1):
    mpr=m2
    y_p1=obs_Pressure2[t]    ##观测值一定要读对，不要读成基质吸力了
    y_obs_p1=pd.DataFrame(y_p1)
    y_obs_p1=y_obs_p1.T
    y_obs_p1=y_obs_p1.values
    y_obs_p1_ave=np.mean(y_obs_p1,axis=1)
    y_obs_prediction=np.zeros_like(y_obs_p1)
    for i in xrange(N):
        y_obs_p1_error[:,i]=y_obs_p1[:,i]-y_obs_p1_ave[:]
        kl_error[:,i]=mpr[:,i]-parY[:,t-2]
    Cd=np.dot(y_obs_p1_error,y_obs_p1_error.T)/(N-1)
    Cm=np.dot(kl_error,kl_error.T)/(N-1)
    m1=mpr
    Imax=1
    
    while True:
        print '{0}_time_{1}_iteration'.format(t,Imax)
        mat_to_arr=lambda x:np.array(x)
        m1,y_obs_prediction,y_obs_p1=map(mat_to_arr,[m1,y_obs_prediction,y_obs_p1])
        m1_average=np.mean(m1,axis=1)
        for i in xrange(N):
            m1_error[:,i]=m1[:,i]-m1_average[:]
            
        ks=np.exp(np.log(ki_mean)+np.dot(Fi[:,:kl_term],m1))
        
        for i in xrange(N):
            time_modify(t,i,filename_time,root_directory)
            ki_modify=ks[:,i]
            ki_modify=para_keys_modify(ki_modify)
            write_ki(Nod_num,ki_modify,i,filename_KI,root_directory)  #para为初始参数，样本之间存在随机扰动

        jobs=[]
        for i in xrange(N):
            p=multiprocessing.Process(target=runexe,args=(i,root_directory,))
            jobs.append(p)
            p.start()
        for j in jobs:
            j.join()       
        
        y_obs_prediction=np.array(y_obs_prediction)
        for k in xrange(N):
            Obs_p1=read_obs(Nod_num,obs_Num,k,root_directory,filename_result)
            for ii in xrange(obs_num):
                y_obs_prediction[ii][k]=Obs_p1[ii]  
        y_obs_prediction_1=y_obs_prediction
         
       
        y_obs_prediction_average=np.mean(y_obs_prediction,axis=1)
        for i in xrange(N):
            y_obs_prediction_error[:,i]=y_obs_prediction[:,i]-y_obs_prediction_average[:]
        
        wn_M=np.linalg.pinv(m1_error)
        G=np.dot(y_obs_prediction_error,wn_M)
        arr_to_mar=lambda x: np.matrix(x)
        mpr,m1,Cm,G,Cd,y_obs_prediction,y_obs_p1=map(arr_to_mar,[mpr,m1,Cm,G,Cd,y_obs_prediction,y_obs_p1])
        m2=beta*mpr+(1-beta)*m1-beta*Cm*(G.T)*(np.linalg.inv(Cd+G*Cm*G.T))*(y_obs_prediction-y_obs_p1-G*(m1-mpr))
        S1=np.trace(((y_obs_prediction-y_obs_p1).T)*np.linalg.inv(Cd)*(y_obs_prediction-y_obs_p1))
        m2=np.array(m2)
        y_obs_prediction=np.array(y_obs_prediction)
        
        
        ks2=np.exp(np.log(ki_mean)+np.dot(Fi[:,:kl_term],m2))
        for i in xrange(N):
            ki_modify=ks2[:,i]
            ki_modify=para_keys_modify(ki_modify)
            write_ki(Nod_num,ki_modify,i,filename_KI,root_directory)  #para为初始参数，样本之间存在随机扰动

        
        jobs=[]
        for i in xrange(N):
            p=multiprocessing.Process(target=runexe,args=(i,root_directory,))
            jobs.append(p)
            p.start()
        for j in jobs:
            j.join()
        
        for k in xrange(N):
            Obs_p1=read_obs(Nod_num,obs_Num,k,root_directory,filename_result)
            for ii in xrange(obs_num):
                y_obs_prediction[ii][k]=Obs_p1[ii]  #组成观测点处的预测值矩阵
        # np.savetxt('y_obs_prediction_{0}_2.txt'.format(t),y_obs_prediction)
        y_obs_prediction=np.matrix(y_obs_prediction)
        
        S2=np.trace(((y_obs_prediction-y_obs_p1).T)*np.linalg.inv(Cd)*(y_obs_prediction-y_obs_p1))
        print 'S1:',S1
        print 'S2:',S2
        print 'max(m2-m1):',np.max(np.abs(m2-m1))
        print 'S2-S1:',S2-S1
        print 'eps2*S1',eps2*S1,'\n'
        
        if S2<S1:
            if (np.max(np.abs(m2-m1))<eps1 or S1-S2<eps2*S1 or Imax>3):
                break
            m1=m2
            beta=2*beta
            if beta>1:
                beta=1
        else:
            beta=0.5*beta
        Imax+=1
        if (Imax>5):
            break
    iter_num.append(Imax*np.float(t)/time_step)
    # np.savetxt('t_{0}.txt'.format(t),m2)
    parY[:,t-1]=np.average(m2,axis=1)
    updated_para_each_time=np.log(ki_mean)+np.dot(Fi[:,:kl_term],m2)
    updated_para_average=np.mean(updated_para_each_time,axis=1)
    updated_para[:,t-1]=updated_para_average
    spread[t-1]=np.sqrt(np.average(np.var(updated_para_each_time,axis=1)))
    rmse[t-1]=np.sqrt(np.sum((updated_para_average-np.squeeze(Y_true))**2)/(Nod_num-1))

np.savetxt('updated_para_{0}_{1}.txt'.format(N,index),updated_para[:,-1])
np.savetxt('updated_para_sample_{0}_{1}.txt'.format(N,index),updated_para_each_time)
iter_num=np.array(iter_num) 
total_iter_cost=np.sum(iter_num*2*N)
print 'average iterion cost',total_iter_cost/time_step
np.savetxt('average_iter_cost_{0}_{1}.txt'.format(N,index),np.array([total_iter_cost/time_step]))
# fig,ax=plt.subplots(1,2,figsize=(12,5))
# ax[0].plot(range(len(rmse)),rmse)
# ax[1].plot(range(len(spread)),spread)
# ax[0].set_xlabel('time')
# ax[0].set_ylabel('RMSE')
# ax[0].set_ylim(0,np.max(rmse)+0.1)
# ax[1].set_xlabel('time')
# ax[1].set_ylabel('Spread')
# ax[1].set_ylim(0,np.max(spread)+0.1)
# fig,ax=plt.subplots()
# ax.scatter(range(1,len(iter_num)+1),iter_num*2*N,s=40,c='r',marker='o')
# print 'average iterations:',iter_num*2*N
# ax.plot(range(1,len(iter_num)+1),iter_num*2*N)
# ax.set_ylabel('Number of model evaluations')
# ax.set_xlabel('time')
# # ax.set_ylim(0,7*2*N)
# ax.set_xlim(0,21)
# plt.show()