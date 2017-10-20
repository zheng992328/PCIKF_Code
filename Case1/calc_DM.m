load para_true.mat
load updated_para_PCIKF.mat
time_step=10;
total_index=10;
obs=ogs_call(para_true,999,time_step);
prediction_PCIKF=ogs_call(updated_para_PCIKF,999,time_step);
DM_PCIKF=sqrt(mean((obs-prediction_PCIKF).^2))

DM_EnRML_80_mat=zeros(total_index);
for i=1:total_index
    load( ['updated_para_EnRML_',num2str(i),'_',num2str(80),'.mat'])
    prediction_EnRML=ogs_call(updated_para_EnRML,999,time_step);
    DM_EnRML_80_mat=sqrt(mean((obs-prediction_EnRML).^2));
end
DM_EnRML_80=mean(DM_EnRML_80_mat)

DM_EnRML_100_mat=zeros(total_index);
for i=1:total_index
    load( ['updated_para_EnRML_',num2str(i),'_',num2str(100),'.mat'])
    prediction_EnRML=ogs_call(updated_para_EnRML,999,time_step);
    DM_EnRML_100_mat=sqrt(mean((obs-prediction_EnRML).^2));
end
DM_EnRML_100=mean(DM_EnRML_100_mat)    