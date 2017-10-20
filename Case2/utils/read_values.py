#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'zq'

import pandas as pd
import numpy as np
import os



##定义函数，读取观测点处对应的预测值
def read_obs(Nod_num,obs_Num,i,root_directory,filename_result):
    '''
    读取当前时刻观测点处对应的预测值
    '''
    gas_file='gas_{0}'.format(i)
    args_obs=os.path.join(root_directory,gas_file,filename_result)
    with open(args_obs,'r') as f:
        Nod_inf=f.readlines()
    last=Nod_inf[-(1600+Nod_num):-1600]  #1600是elements的个数
    for i in xrange(Nod_num):
        last[i]=last[i].split()
    last=pd.DataFrame(last)
    obs_prediction=last.ix[obs_Num,[4]]
    obs_prediction=obs_prediction.values
    obs_prediction=np.float64(obs_prediction)
    obs_p1=obs_prediction[:,0]
    return obs_p1

def read_obs_alltime(t,i,obs_Num,Nod_num,root_directory,filename_result):
    '''
    t表示取出1~t时刻的所有观测点的预测值
    '''
    gas_file='gas_{0}'.format(i)
    args_domain=os.path.join(root_directory,gas_file,filename_result)
    with open(args_domain,'r') as f:
        content=f.readlines()
    x=np.zeros((0,5))   #准备提取tec文件中每个时间步输出的前5列数据（x,y,z,p1,p2）
    for i in xrange(t+1):
        content1=content[3284*i:3284*i+Nod_num+3]  #跳过第0个时间步，第一个时间步从第3284行开始，854会随网格划分格点及单元个数而改变
        content1=content1[3:]              #踢掉tec前面3行头数据
        for i in xrange(len(content1)):
            content1[i]=content1[i].split()
        content1=pd.DataFrame(content1)
        content_1=content1.ix[obs_Num,[0,1,2,3,4]]
        value=content_1.values
        x=np.vstack((x,value))

    x=x[len(obs_Num):]   ##把0时刻的踢掉
    x=np.array(x)
    x=np.float64(x)
    return x[:,4]
###读取所有的p2值
def read_p(Nod_num,i,filename_result,root_directory):
   gas_file='gas_{0}'.format(i)
   p_result=os.path.join(root_directory,gas_file,filename_result)
   with open(p_result,'r') as f:
       Nod_inf=f.readlines()
   last=Nod_inf[-(1600+Nod_num):-1600]
   head=[None]*Nod_num
   for i in xrange(len(last)):
       line=last[i].split()
       head[i]=line[4:5]

   head=np.array(head)
   head=np.float64(head)
   return head