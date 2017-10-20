# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
from write_values import write_ki
import numpy as np
import win32api
from KLexpansion import KLexpansion
import time,pp
import os
##提取一个点（0.75,0.4）处的20个时间步的观测值
current_directory=os.getcwd()
root_directory=os.path.dirname(current_directory)
path_for_add='true_obs\obs_20.txt'
point_args=os.path.join(root_directory,path_for_add)
obs_values=pd.read_csv(point_args,sep=' ',names=[i for i in range(5)])
line=[5+16*i for i in range(20)]
values_for_map=obs_values.ix[line,[4]]
values_for_map=values_for_map.values
filename_KI='water_gas_KI.direct'
ogs_para='water_gas'
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
kl_term=20
Fi=KLexpansion(m,n,x,y,deltax,deltay,sigma)
# ##一个文件夹提取出一根线(0.75,0.4)的20个时间步的变化图）
def get_obs_updated(i):
    line_num_updated=[3918+3284*ii for ii in range(20)]
    filename='gas_{0}\water_gas_domain_quad.tec'.format(i)
    args_result=os.path.join(root_directory,filename)
    with open(args_result,'r') as f:
        content=f.readlines()
    for i in xrange(len(content)):
        content[i]=content[i].split()
    content=np.array(content)
    ff=lambda x:np.float64(x)
    obs_updated=content[line_num_updated]
    obs_updated=map(ff,obs_updated)
    obs_up=np.zeros((20,5))
    for m in range(20):
        for n in range(5):
            obs_up[m,n]=obs_updated[m][n]
    result=obs_up[:,4]
    return result


N=100
Nod_num=1681
args_para_updated=os.path.join(root_directory,'t_20.txt')
args_para_initial=os.path.join(root_directory,'para_initial.txt')
kl_updated=np.loadtxt(args_para_updated)
para_updated=np.log(ki_mean)+np.dot(Fi[:,:kl_term],kl_updated)
para_initial=np.loadtxt(args_para_initial)
##将初始参数带进去，计算得到结果
def runexe(i,root_directory):
    path_1=root_directory
    path_2='gas_{0}\ogs.exe'.format(i)
    path_3='gas_{0}'.format(i)
    args_exe=os.path.join(path_1,path_2)
    args=os.path.join(path_1,path_3)
    win32api.ShellExecute(0,'open',args_exe,'water_gas',args,0)



  
for i in xrange(N):
    para_init=para_initial[:,i]
    para_initial_=np.exp(para_init)
    write_ki(1681,para_initial_,i,filename_KI,root_directory)

ppservers=()
job_server=pp.Server(4,ppservers=ppservers)
Ne=tuple(range(N))
jobs=[(i,job_server.submit(runexe,(i,root_directory),(),('win32api',))) for i in Ne]
for i,job in jobs:
    job()

while True:
    back_process=os.popen('tasklist |findstr /i ogs.exe').readlines()
    if back_process==[]:
        break
os.popen('taskkill /f /fi "imagename eq python.exe" /fi "memusage lt 100000"')

values_initial=np.zeros((20,N))
for i in xrange(N):
    obs_initial=get_obs_updated(i)
    for m in xrange(len(obs_initial)):
        values_initial[m,i]=obs_initial[m]

#
###将更新后的参数带进去，计算得到结果

for i in xrange(N):
    para_up=para_updated[:,i]
    para_updated_=np.exp(para_up)
    write_ki(1681,para_updated_,i,filename_KI,root_directory)
    
ppservers=()
job_server=pp.Server(4,ppservers=ppservers)
Ne=tuple(range(N))
jobs=[(i,job_server.submit(runexe,(i,root_directory),(),('win32api',))) for i in Ne]
for i,job in jobs:
    job()


while True:
    back_process=os.popen('tasklist |findstr /i ogs.exe').readlines()
    if back_process==[]:
        break

os.popen('taskkill /f /fi "imagename eq python.exe" /fi "memusage lt 100000"')


values_updated=np.zeros((20,N))
for i in xrange(N):
    obs_update=get_obs_updated(i)
    for m in xrange(len(obs_update)):
        values_updated[m,i]=obs_update[m]


x=[i+1 for i in range(20)]
N=100
#plt.figure(figsize=(6,5))
fig,ax=plt.subplots(1,2,figsize=(13,6))
for k in range(N+1):
    if 1<k<=N-1:
        ax[0].plot(x,values_initial[:,k-1],'b',alpha=0.4)     
        ax[1].plot(x,values_updated[:,k-1],'y',alpha=0.3)
    if k==N:
        ax[0].plot(x,values_initial[:,k-1],'b',alpha=0.4,label='initial')    
        ax[1].plot(x,values_updated[:,k-1],'y',alpha=0.8,label='updated')
    if k==0:
        ax[0].scatter(x,values_for_map,color='r',marker='x',s=60,label='observation')
        ax[1].scatter(x,values_for_map,color='r',marker='x',s=60,label='observation')
      
ax[0].set_xlim(0,20)
ax[1].set_xlim(0,20)
ax[0].set_xlabel('time step')
ax[1].set_xlabel('time step')
ax[0].set_ylabel('Pressure(Pa)')
ax[0].legend(loc='best')
ax[1].legend(loc='best')
plt.show()





