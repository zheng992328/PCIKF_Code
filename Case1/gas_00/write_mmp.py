# -*- coding: utf-8 -*-



import os
import numpy as np
import scipy.io as si
 
def write_mmp(para,filename):
    # gas_file='gas_{0}'.format(i)
    # path=os.path.join(root_directory,gas_file,filename)
    para=np.squeeze(para)
    with open(filename,'r') as f:
        content=f.readlines()
    content_for_modify_1=content[9].split()  ##ä¿®æ¹ç¬¬å è¡æ°é?¦çµæ´»è°æ´
    content_for_modify_2=content[27].split()
    content_for_modify_3=content[45].split()
    # content_for_modify_4=content[63].split()
    content_for_modify_1[-1]=str(para[0])
    content_for_modify_2[-1]=str(para[1])
    content_for_modify_3[-1]=str(para[2])
    # content_for_modify_4[-1]=str(para[3])
    content_for_modify_1=' '.join(element for element in content_for_modify_1)+'\n'
    content_for_modify_2=' '.join(element for element in content_for_modify_2)+'\n'
    content_for_modify_3=' '.join(element for element in content_for_modify_3)+'\n'
    # content_for_modify_4=' '.join(element for element in content_for_modify_4)+'\n'
    content[9]=content_for_modify_1
    content[27]=content_for_modify_2
    content[45]=content_for_modify_3
    # content[63]=content_for_modify_4
    with open(filename,'w') as f:
        for line in content:
            f.write(line)
  
if __name__=='__main__':
    current_directory=os.getcwd()
    root_directory=os.path.dirname(current_directory)
    filename='first_case.mmp'
#    para=np.array([1e-12,2e-12,3e-12])
    para=si.loadmat('x.mat')
    para=para['x'].astype(float)
    para=np.exp(para)
    write_mmp(para,filename)
    