# -*- coding: utf-8 -*-
# import win32api
import os
import subprocess
# def runexe(i,root_directory):
    # path_1=root_directory
    # path_2='gas_{0}\ogs.exe'.format(i)
    # path_3='gas_{0}'.format(i)
    # args_exe=os.path.join(path_1,path_2)
    # args=os.path.join(path_1,path_3)
    # win32api.ShellExecute(0,'open',args_exe,'water_gas',args,0)
    
def runexe(i,root_directory):
    path_1=root_directory
    path_2='gas_{0}/ogs'.format(i)
    path_3='gas_{0}/water_gas'.format(i)
    args_exe=os.path.join(path_1,path_2)
    para=os.path.join(path_1,path_3)
    pp=subprocess.Popen((args_exe,para),stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    pp.communicate()
    pp.wait()