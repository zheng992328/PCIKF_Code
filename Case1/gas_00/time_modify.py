#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'zq'

import sys
import numpy as np
def time_modify(t,filename_time):
    t=np.int(t)
    with open(filename_time,'r') as f:
        content=f.readlines()
        line_modify_1=content[5]
        line_modify_2=content[7]
        modify_1=line_modify_1.split()
        modify_1[0]='{0}'.format(t)
        line_modify_1=' '.join(modify_1[i] for i in range(len(modify_1)))+'\n'
        modify_2=line_modify_2.split()
        modify_2[0]='{0}'.format(t*600)
        line_modify_2=' '.join(modify_2[i] for i in range(len(modify_2)))+'\n'
        content[5]=line_modify_1
        content[7]=line_modify_2

    with open(filename_time,'w') as f:
        for line in content:
            f.write(line)

if __name__=='__main__':
    t=sys.argv[1] #这里的t进来是字符串，所以函数当中要对其进行Int转换
    time_modify(t,'first_case.tim')

