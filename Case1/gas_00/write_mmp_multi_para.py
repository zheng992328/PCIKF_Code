# -*- coding: utf-8 -*-



import os
import numpy as np
 
def write_mmp(para,filename):
   
    with open(filename,'r') as f:
        content=f.readlines()
    content_for_modify_1=content[9].split()  ##ä¿®æ¹ç¬¬å è¡æ°é?¦çµæ´»è°æ´
    content_for_modify_2=content[14].split()  ##ä¿®æ¹ç¬¬å è¡æ°é?¦çµæ´»è°æ´
    content_for_modify_3=content[15].split()
    content_for_modify_4=content[17].split()
    
    content_for_modify_5=content[27].split()
    content_for_modify_6=content[32].split()  ##ä¿®æ¹ç¬¬å è¡æ°é?¦çµæ´»è°æ´
    content_for_modify_7=content[33].split()
    content_for_modify_8=content[35].split()
    
    content_for_modify_9=content[45].split()
    content_for_modify_10=content[50].split()  ##ä¿®æ¹ç¬¬å è¡æ°é?¦çµæ´»è°æ´
    content_for_modify_11=content[51].split()
    content_for_modify_12=content[53].split()
    
    content_for_modify_1[-1]=str(para[0])
    content_for_modify_2[1]=str(para[1])
    content_for_modify_2[2]=str(para[2])
    content_for_modify_2[3]=str(para[3])
    content_for_modify_3[1]=str(1-para[2])
    content_for_modify_3[2]=str(1-para[1])
    content_for_modify_3[3]=str(para[3])
    content_for_modify_4[1]=str(para[4])
    content_for_modify_4[2]=str(para[1])
    content_for_modify_4[3]=str(para[2])
    content_for_modify_4[4]=str(para[3])
    

    content_for_modify_5[-1]=str(para[5])
    content_for_modify_6[1]=str(para[6])
    content_for_modify_6[2]=str(para[7])
    content_for_modify_6[3]=str(para[8])
    content_for_modify_7[1]=str(1-para[7])
    content_for_modify_7[2]=str(1-para[6])
    content_for_modify_7[3]=str(para[8])
    content_for_modify_8[1]=str(para[9])
    content_for_modify_8[2]=str(para[6])
    content_for_modify_8[3]=str(para[7])
    content_for_modify_8[4]=str(para[8])
    

    content_for_modify_9[-1]=str(para[10])
    content_for_modify_10[1]=str(para[11])
    content_for_modify_10[2]=str(para[12])
    content_for_modify_10[3]=str(para[13])
    content_for_modify_11[1]=str(1-para[12])
    content_for_modify_11[2]=str(1-para[11])
    content_for_modify_11[3]=str(para[13])
    content_for_modify_12[1]=str(para[14])
    content_for_modify_12[2]=str(para[11])
    content_for_modify_12[3]=str(para[12])
    content_for_modify_12[4]=str(para[13])
    
    
    for i in range(1,13):
        variable='content_for_modify_{0}'.format(i)
        exec(variable+"=' '.join(element for element in "+ variable +")+'\\n'")
    
    
    content[9]=content_for_modify_1
    content[14]=content_for_modify_2
    content[15]=content_for_modify_3
    content[17]=content_for_modify_4
    
    content[27]=content_for_modify_5
    content[32]=content_for_modify_6
    content[33]=content_for_modify_7
    content[35]=content_for_modify_8
    
    content[45]=content_for_modify_9
    content[50]=content_for_modify_10
    content[51]=content_for_modify_11
    content[53]=content_for_modify_12
    
    with open(filename,'w') as f:
        for line in content:
            f.write(line)
  
if __name__=='__main__':
    import scipy.io as si
    filename='first_case.mmp'
    para=si.loadmat('x.mat')
    para=para['x'].astype(float)
    para=np.squeeze(para)
    para[0]=np.exp(para[0])
    para[5]=np.exp(para[5])
    para[10]=np.exp(para[10])
    para[4]=1000*para[4]
    para[9]=1000*para[9]
    para[14]=1000*para[14]
    write_mmp(para,filename)
    
    
    
   