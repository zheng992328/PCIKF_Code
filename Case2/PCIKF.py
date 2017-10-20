# -*- coding: utf-8 -*-

import numpy as np
import os
import scipy as sp
import matplotlib.pyplot as plt
import time
from utils.KLexpansion import KLexpansion
from utils.write_values import write_ki
from utils.runexe import runexe
from utils.para_key_mod import para_keys_modify
from utils.read_values import read_obs_alltime
from utils.get_PCE_coef_parallel import get_PCE_coef
from utils.time_modify import time_modify

start=time.time()
Obs_Num=16
R=np.eye(Obs_Num)
varR=80
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
m=int(y/dy+1)
n=int(x/dx+1)
para_for_ogs='water_gas'
filename_result='water_gas_domain_quad.tec'
filename_KI='water_gas_Ki.direct'
filename_time='water_gas.tim'
root_directory=os.getcwd()
Fi=KLexpansion(m,n,x,y,deltax,deltay,sigma)
beta=0.05
eps1=0.03
eps2=1e-2

f=open('log.txt','w')

##产生观测值
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
y_obs=np.zeros((time_step*Obs_Num,1+2*kl_term))
y_obs[:,0]=obs

C_P0=Y_mean
C_P1=np.zeros((Nod_num,2,kl_term))
C_P1[:,0,:]=Fi[:,:kl_term]
# C_P2=np.zeros((Nod_num,2,2,kl_term,kl_term))
C_P=np.hstack((C_P0,C_P1.reshape((Nod_num,2*kl_term))))  #参数先验
Num=np.zeros(time_step)
parY=np.zeros((Nod_num,time_step))
rmse=[]
spread=[]
rs=[]
iter_num=[]
time_sum=np.sum(np.arange(1,21))
total_Imax=[]  #统计每个时间步迭代次数，以便计算平均调用次数，与样本个数等价
for t in xrange(1,time_step+1):
    print 'the {0} time_step'.format(t)
    f.write('the {0} time_step'.format(t)+'\n')
    
    C_Prior=C_P  #每一步更新出来用C_P接，并且将其作为下一步的先验，同时也是该步迭代的初始值
    Cm=np.dot(C_Prior[:,1:],C_Prior[:,1:].T)
    Cd=np.eye(Obs_Num)*varR
    C_P_1=C_Prior
          
        
    Imax=1
    while True:
        print 'this is {0}_time {1}_iteration'.format(t,Imax)
        f.write('this is {0}_time {1}_iteration'.format(t,Imax)+'\n')
        # np.savetxt('{0}_time_{1}_iter_C_P.txt'.format(t,Imax),C_P_1)
        C_O0,C_O1=get_PCE_coef(kl_term,Nod_num,C_P0,C_P1,obs_Num,sigma,ki_mean,deltax,deltay,dx,dy,x,y,para_for_ogs,filename_result,filename_KI,filename_time,root_directory,out_node,t)
        C_O=np.hstack((C_O0.reshape((-1,1)),C_O1.reshape((Obs_Num,2*kl_term))))
        G=np.dot(C_O[:,1:],np.linalg.pinv(C_P_1[:,1:]))
        [C_P,C_Prior,C_P_1,Cm,G,Cd,C_O,y_obs]=map(lambda x:np.matrix(x),[C_P,C_Prior,C_P_1,Cm,G,Cd,C_O,y_obs])
        
        C_P=beta*C_Prior+(1-beta)*C_P_1-beta*Cm*(G.T)*(np.linalg.inv(Cd+G*Cm*G.T))*(C_O-y_obs[Obs_Num*(t-1):Obs_Num*t,:]-G*(C_P_1-C_Prior))
        S1=np.trace(((C_O-y_obs[Obs_Num*(t-1):Obs_Num*t,:]).T)*np.linalg.inv(Cd)*(C_O-y_obs[Obs_Num*(t-1):Obs_Num*t,:]))
        [C_P,C_O]=map(lambda x:np.array(x),[C_P,C_O])
        # np.savetxt('{0}_time_{1}_iter_C_O.txt'.format(t,Imax),C_O)
        
        ##将每一次迭代之后计算得到的即将继续迭代的值保存下来
        # np.savetxt('{0}_time_{1}_iter_C_P.txt'.format(t,Imax),C_P[:,0])
        
        # 更新之后的计算
        C_P0=C_P[:,0]
        C_P1=C_P[:,1:1+2*kl_term].reshape((Nod_num,2,kl_term))
        # C_P2=C_P[:,1+2*kl_term:].reshape((Nod_num,2,2,kl_term,kl_term))
        C_O0,C_O1=get_PCE_coef(kl_term,Nod_num,C_P0,C_P1,obs_Num,sigma,ki_mean,deltax,deltay,dx,dy,x,y,para_for_ogs,filename_result,filename_KI,filename_time,root_directory,out_node,t)
        C_O=np.hstack((C_O0.reshape((-1,1)),C_O1.reshape((Obs_Num,2*kl_term))))
        C_O=np.matrix(C_O)

        S2=np.trace(((C_O-y_obs[Obs_Num*(t-1):Obs_Num*t,:]).T)*np.linalg.inv(Cd)*(C_O-y_obs[Obs_Num*(t-1):Obs_Num*t,:]))

        
        
        print 'S1:',S1
        print 'S2:',S2
        print 'max(C_P-C_P_1):',np.max(np.abs(C_P-C_P_1))
        print 'S2-S1:',S2-S1
        print 'eps2*S1',eps2*S1,'\n'
        f.write('S1:{0}'.format(S1)+'\n')
        f.write('S2:{0}'.format(S2)+'\n')
        f.write('max(C_P-C_P_1):{0}'.format(np.max(np.abs(C_P-C_Prior)))+'\n')
        f.write('S2-S1:{0}'.format(S2-S1)+'\n')
        f.write('eps2*S1:{0}'.format(eps2*S1)+'\n')
        
        if S2<S1:
            if (np.max(np.abs(C_P-C_P_1))<eps1 or S1-S2<eps2*S1 or Imax>3):
                f.write('\n')
                break
            C_P_1=C_P
            C_P0=C_P_1[:,0]
            C_P1=C_P_1[:,1:1+2*kl_term].reshape((Nod_num,2,kl_term))
            # C_P2=C_P_1[:,1+2*kl_term:].reshape((Nod_num,2,2,kl_term,kl_term))
            beta=2*beta
            f.write('beta:{0}'.format(beta))
            f.write('\n')
            if beta>1:
                beta=1
        else:
            beta=0.5*beta   # 如果S2一直大于S1，会一直走else走下去吗，好像这个逻辑也有问题
            f.write('beta:{0}'.format(beta))
            f.write('\n')
            
        Imax+=1
        if (Imax>5):
            break
    total_Imax.append(Imax)
    iter_num.append(Imax*np.float(t)/time_step)
    parY[:,t-1]=C_P0
    
    rmse_=np.sqrt(np.sum((C_P[:,0]-np.squeeze(Y_true))**2)/(Nod_num-1))
    rmse.append(rmse_)
    temp_1=C_P1.reshape((Nod_num,2*kl_term))
    # temp_2=C_P2.reshape((Nod_num,(2*kl_term)**2))
    # temp_3=np.hstack((temp_1,temp_2))
    spread_=np.sqrt(np.average(np.sum(temp_1**2,axis=1)))
    spread.append(spread_)
    rs_=2*np.sum(C_P[:,0]*(C_P[:,0]-np.squeeze(Y_true)))/Nod_num+np.sum(np.squeeze(Y_true)**2-C_P[:,0]**2)/Nod_num
    rs.append(rs_)
iter_num=np.array(iter_num)
np.savetxt('total_Imax.txt',total_Imax)
np.savetxt('pcrml_rmse.txt',rmse)
np.savetxt('pcrml_spread.txt',spread)
np.savetxt('rmse_spread.txt',rs)
total_iter_cost=np.sum(iter_num*2*(1+2*kl_term))
average_iter_cost=total_iter_cost/time_step
np.savetxt('average_iter_cost.txt',np.array([average_iter_cost]))

updated_para=parY[:,-1]
np.savetxt('updated_para.txt',updated_para)
f.close()

fig,ax=plt.subplots(1,2,figsize=(12,5))
ax[0].plot(range(len(rmse)),rmse)
ax[1].plot(range(len(spread)),spread)
ax[0].set_xlabel('time')
ax[0].set_ylabel('RMSE')
ax[0].set_ylim(0,np.max(rmse)+0.1)
ax[1].set_xlabel('time')
ax[1].set_ylabel('Spread')
ax[1].set_ylim(0,np.max(spread)+0.1)
fig.savefig('RS.jpeg')
fig,ax=plt.subplots()
ax.scatter(range(1,len(iter_num)+1),iter_num*2*(1+2*kl_term),s=40,c='r',marker='o')
ax.plot(range(1,len(iter_num)+1),iter_num*2*(1+2*kl_term))
ax.set_ylabel('Number of model evaluations')
ax.set_xlabel('time')
# ax.set_ylim(0,2*(1+2*kl_term))
ax.set_xlim(0,21)
fig.savefig('iteration.jpeg')
plt.show()



