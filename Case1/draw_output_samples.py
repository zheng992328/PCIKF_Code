# -*- coding: utf-8 -*-


#为了保持与第二个案例的绘图一致性，用matlab版本的draw_output_samples.m计算好结果，再用python绘图



import numpy as np
import matplotlib.pyplot as plt

N=100
time_index=10
time_step=10
obs_num=6

y_obs_prediction=np.loadtxt('y_obs_ensem_{0}_time_{1}.txt'.format(N,time_index))
obs_mean=np.loadtxt('y_obs_mean_{0}_time_{1}.txt'.format(N,time_index))
obs_true=np.loadtxt('true_obs_{0}_time_{1}.txt'.format(N,time_index))


#提取3个观测点的，每个时间步的数据
obs_index_ensem=[1,3,5]
fig,ax=plt.subplots(1,3,figsize=(19,6),squeeze=True)
for obs_index in range(len(obs_index_ensem)):
    obs_ensemble=np.zeros((time_step,N))
    for i in xrange(N):
        for t_index in xrange(time_step):
            obs_ensemble[t_index,i]=y_obs_prediction[obs_index_ensem[obs_index]+t_index*obs_num-1,i]
            
    obs_mean_selected=np.zeros(time_step)      
    for t_index in xrange(time_step):
        obs_mean_selected[t_index]=obs_mean[obs_index_ensem[obs_index]+t_index*obs_num-1]
        
    obs_true_selected=np.zeros(time_step)      
    for t_index in xrange(time_step):
        obs_true_selected[t_index]=obs_true[obs_index_ensem[obs_index]+t_index*obs_num-1]
    
    for k in xrange(N):
        if k==N-1:
            ax[obs_index].plot(range(time_step),obs_ensemble[:,k],c='b',alpha=0.2,label='Realizations')
            ax[obs_index].plot(range(time_step),obs_mean_selected,c='r',lw=1.5,alpha=1,label='Mean estimation')
            ax[obs_index].scatter(range(time_step),obs_true_selected,c='k',s=22,label='Observations',zorder=10)
        ax[obs_index].plot(range(time_step),obs_ensemble[:,k],c='b',alpha=0.1)
        min_y=np.int32(min(ax[obs_index].get_ylim()))
        max_y=np.int32(max(ax[obs_index].get_ylim()))
        ax[obs_index].set_ylim(min_y,max_y)
        ax[obs_index].set_yticks(np.linspace(min_y,max_y,6))
        ax[obs_index].set_yticklabels(np.int32(np.linspace(min_y,max_y,6)))
        
        
    ax[obs_index].plot(range(time_step),obs_mean_selected,c='r',lw=1.5,alpha=1)
    ax[obs_index].scatter(range(time_step),obs_true_selected,c='k',s=22)
    ax[obs_index].set_xlabel('Time step',fontsize=15)
    ax[obs_index].set_ylabel('Pressure(Pa)',fontsize=15)
    ax[obs_index].set_xlim(-1,11)
    plt.setp(ax[obs_index].get_xticklabels(),fontsize=14)
    plt.setp(ax[obs_index].get_yticklabels(),fontsize=14)




    
# ax[2].set_ylim(100800,101400)
# ax[2].set_yticks(np.linspace(100800,101400,6))
# ax[2].set_yticklabels(np.int32(np.linspace(100800,101400,6)))
plt.legend(loc='best')
plt.tight_layout() 
fig.subplots_adjust(wspace=0.25)
fig.savefig('PCIKF_fitting_map_time_{0}.jpeg'.format(time_index))  
plt.show()







