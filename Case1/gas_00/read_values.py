#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'zq'

import pandas as pd
import numpy as np
import os,sys,re



##定义函数，读取观测点处对应的预测值
def read_obs(Nod_num,obs_Num,i,root_directory,filename_result):
    '''
    读取当前时刻观测点处对应的预测值
    '''
    # gas_file='gas_{0}'.format(i)
    # args_obs=os.path.join(root_directory,gas_file,filename_result)
    with open(filename_result,'r') as f:
        Nod_inf=f.readlines()
    last=Nod_inf[-(900+Nod_num):-900]  #900是elements的个数
    for i in xrange(Nod_num):
        last[i]=last[i].split()
    last=pd.DataFrame(last)
    obs_prediction=last.ix[obs_Num,[3]]
    obs_prediction=obs_prediction.values
    obs_prediction=np.float64(obs_prediction)
    obs_p1=obs_prediction[:,0]
    return obs_p1

def read_obs_alltime(t,obs_Num,Nod_num,root_directory,filename_result):
    '''
    t表示取出t-最后时刻的所有观测点的预测值
    '''
    # gas_file='gas_{0}'.format(i)
    # args_domain=os.path.join(root_directory,gas_file,filename_result)

    with open(filename_result,'r') as f:
        content=f.readlines()
    x=np.zeros((0,4))   #准备提取tec文件中每个时间步输出的前5列数据（x,y,z,p1,p2）
    for i in xrange(1,t+1):
        content1=content[1864*i:1864*i+Nod_num+3]  #跳过第0个时间步，第一个时间步从第3284行开始，854会随网格划分格点及单元个数而改变
        content1=content1[3:]              #踢掉tec前面3行头数据
        for i in xrange(len(content1)):
            content1[i]=content1[i].split()
        content1=pd.DataFrame(content1)
        content_1=content1.ix[obs_Num,[0,1,2,3]]
        value=content_1.values
#        print value
        x=np.vstack((x,value))

    x=np.array(x)
    x=np.float64(x)
    return x[:,-1]

   
   
if __name__=='__main__':
    obs_Num=[315,325,335,625,635,645]
    Nod_num=961
    root_directory=os.getcwd()
    end=root_directory[-3:]  ##获得当前路径是第几个gas文件夹
    n=re.findall(r'\d',end)
    i=''.join(element for element in n)
    # root_directory=os.path.dirname(current_directory)
    filename_result='first_case_domain_quad.tec'
#    aa=read_obs_alltime(t,obs_Num,Nod_num,root_directory,filename_result)
    aa=read_obs(Nod_num,obs_Num,i,root_directory,filename_result)
    np.savetxt('data_all_time_{0}.txt'.format(i),aa)
    
    
    
    