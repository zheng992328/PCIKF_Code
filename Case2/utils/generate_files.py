#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'zq'

def remove_generate_files(N):
    import os
    import shutil
    current_directory=os.getcwd()
    root_directory=os.path.dirname(current_directory)
    for i in xrange(N):
        dir_name='gas_{0}'.format(i)
        file_args=os.path.join(root_directory,dir_name)
        old=os.path.join(root_directory,'gas_00')
        if os.path.exists(file_args):
            shutil.rmtree(file_args)
            shutil.copytree(old,file_args)
        else:
            shutil.copytree(old,file_args)

if __name__=='__main__':
    remove_generate_files(100)
 
