
N=100;
time_index=10;
time_step=10;
obs_num=6;

% para_sample=importdata(['./updated_para_sample_EnRML_1_',num2str(N),'_',num2str(time_index),'.mat']);
% para_mean=importdata(['./updated_para_mean_EnRML_1_',num2str(N),'_',num2str(time_index),'.mat']);

[para_sample,para_mean]=generate_PCIKF_post_para(time_index,N);

%计算样本输出
y_obs_ensem=zeros(obs_num*time_step,N);
parfor i=1:N
    y_obs_ensem(:,i)=ogs_call(para_sample(:,i),i,time_step);
end
save(['y_obs_ensem_',num2str(N),'_time_',num2str(time_index),'.txt'],'-ascii','y_obs_ensem');

% load(['y_obs_ensem_',num2str(N),'_time_',num2str(time_index),'.txt'])

%计算均值输出
y_obs_mean=ogs_call(para_mean,1,time_step);
save(['y_obs_mean_',num2str(N),'_time_',num2str(time_index),'.txt'],'-ascii','y_obs_mean')

%得到真实观测值
load para_true.mat
true_obs=ogs_call(para_true,999,time_step);
save(['true_obs_',num2str(N),'_time_',num2str(time_index),'.txt'],'-ascii','true_obs')



%提取3个观测点的每个时间步的观测数据,为了与案例2保持图例一致，转用draw_output_samples.py作图
obs_index_ensem=[1,3,5];

for obs_index=1:length(obs_index_ensem)
    obs_ensemble=zeros(time_step,N);
    for i=1:N
        for t_index=1:time_step
            obs_ensemble(t_index,i)=y_obs_ensem(obs_index_ensem(obs_index)+(t_index-1)*obs_num,i);
        end
    end
    
    obs_mean_selected=zeros(time_step,1);
    for t_index=1:time_step
        obs_mean_selected(t_index)=y_obs_mean(obs_index_ensem(obs_index)+(t_index-1)*obs_num);
    end
    
    obs_true_selected=zeros(time_step,1);
    for t_index=1:time_step
        obs_true_selected(t_index)=true_obs(obs_index_ensem(obs_index)+(t_index-1)*obs_num);
    end
    
    subplot(1,3,obs_index)
    for k=1:N
        plot(1:time_step,obs_ensemble(:,k),'b')
		hold on
    end
    plot(1:time_step,obs_mean_selected,'r-')
	hold on
    plot(1:time_step,obs_true_selected,'k*')
    
end