# -*- coding: utf-8 -*-

##用来绘制t=4,8,12,16,20时刻的参数场分布图，用来和真实参数场对比
import os
import pandas as pd
Nod_num=1681
from para_distribution_map import para_distribution_map
current_directory=os.getcwd()
root_directory=os.path.dirname(current_directory)
root_directory_true_obs=os.path.join(root_directory,'true_obs')

filename='para_updated.txt'
#filename='para_true.txt'
filepath=os.path.join(root_directory,filename)
content=pd.read_csv(filepath,sep=' ',names=['1'])
content=content.values
para=content[:,0]
para_distribution_map(para,Nod_num,root_directory_true_obs)
