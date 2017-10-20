clc,clear;
varR=80.0;
Nod_num=961;
obs_Num=[315,325,335,625,635,645];
Obs_Num=length(obs_Num);
time_step=10;
out_node=Obs_Num;
kl_term=15;
ki_1_mean=1e-12;
ki_2_mean=3e-12;
ki_3_mean=5e-12;
theta_1_r_w=0.24;
theta_2_r_w=0.22;
theta_3_r_w=0.21;
theta_1_s_w=0.63;
theta_2_s_w=0.64;
theta_3_s_w=0.65;
m_v_1=0.33;
m_v_2=0.34;
m_v_3=0.35;
alpha_v_1=0.88;
alpha_v_2=0.9;
alpha_v_3=0.92;
para_mean=[log(ki_1_mean),theta_1_r_w,theta_1_s_w,m_v_1,alpha_v_1,log(ki_2_mean),theta_2_r_w,theta_2_s_w,m_v_2,alpha_v_2,log(ki_3_mean),theta_3_r_w,theta_3_s_w,m_v_3,alpha_v_3];
Fi=diag(repmat([0.6,0.01,0.03,0.01,0.03],1,3));

para_for_ogs='first_case';
filename_result='first_case_domain_quad.tec';
filename_KI='first_case.mmp';
filename_time='first_case.tim';
root_directory=pwd;

beta=0.05;
eps1=3e-6;
eps2=1e-5;

%% 产生观测值
para_true=para_mean'+Fi*randn(kl_term,1);
save 'para_true.mat' para_true
obs=ogs_call(para_true,999,time_step);
y_obs=zeros(time_step*Obs_Num,1+2*kl_term);
y_obs(:,1)=obs;
disp('get obs over');

%% 反演部分
C_P0=para_mean';
C_P1=zeros(kl_term,2,kl_term);
C_P1(:,1,:)=Fi;
C_P=[C_P0,reshape(C_P1,kl_term,2*kl_term)];

Num=zeros(1,time_step);
ParY=zeros(kl_term,time_step);%放更新之后的参数场
RMSE=zeros(1,time_step);
Spread=zeros(1,time_step);
total_Imax=zeros(1,time_step);
iter_num=zeros(1,time_step);
for t=1:time_step
    fprintf('this is the %d time step\n',t);
    C_Prior=C_P;
    Cm=C_Prior(:,2:end)*(C_Prior(:,2:end)');
    Cd=eye(Obs_Num)*varR;
    C_P_1=C_Prior;
    Imax=1;
    
    while 1
        fprintf('this is the %d time %d iteration\n',t,Imax);
        [C_O0,C_O1]=Get_ANOVA_adaptive(root_directory,kl_term,Nod_num,C_P0,C_P1,t,obs_Num,Obs_Num);
        C_O=[C_O0,reshape(C_O1,Obs_Num,2*kl_term)];
        G=C_O(:,2:end)*pinv(C_P_1(:,2:end));
        C_P=beta*C_Prior+(1-beta)*C_P_1-beta*Cm*(G')*(pinv(Cd+G*Cm*G'))*(C_O-y_obs(Obs_Num*(t-1)+1:Obs_Num*t,:)-G*(C_P_1-C_Prior));
        S1=trace((C_O-y_obs(Obs_Num*(t-1)+1:Obs_Num*t,:))'*pinv(Cd)*(C_O-y_obs(Obs_Num*(t-1)+1:Obs_Num*t,:)));
        
        % 更新之后
        C_P0=C_P(:,1);
        C_P1=reshape(C_P(:,2:1+2*kl_term),kl_term,2,kl_term);
        [C_O0,C_O1]=Get_ANOVA_adaptive(root_directory,kl_term,Nod_num,C_P0,C_P1,t,obs_Num,Obs_Num);
        C_O=[C_O0,reshape(C_O1,Obs_Num,2*kl_term)];
        S2=trace((C_O-y_obs(Obs_Num*(t-1)+1:Obs_Num*t,:))'*pinv(Cd)*(C_O-y_obs(Obs_Num*(t-1)+1:Obs_Num*t,:)));
        
        fprintf('S1: %d\n',S1);
        fprintf('S2: %d\n',S2);
        fprintf('max(C_P-C_P_1): %d\n',max(max(abs(C_P-C_P_1))));
        fprintf('S2-S1: %d\n',S2-S1);
        fprintf('eps2*S1: %d\n',eps2*S1);
        if S2<S1
            if max(max(abs(C_P-C_P_1)))<eps1 | S1-S2<eps2*S1 | Imax>2
                break
            end
            C_P_1=C_P;
            C_P0=C_P_1(:,1);
            C_P1=reshape(C_P_1(:,2:1+2*kl_term),kl_term,2,kl_term);
            beta=2*beta;
            if beta>1
                beta=1;
            end
        else
            beta=0.5*beta;           
        end
        Imax=Imax+1;
        if Imax>3
            break
        end
    end
    total_Imax(t)=Imax;
    iter_num(t)=Imax*t/time_step;
    ParY(:,t)=C_P0;
    RMSE(t)=sqrt(mean((C_P0-para_true).^2));
    temp=reshape(C_P1,kl_term,2*kl_term);
    Spread(t)=sqrt(mean(sum(temp.^2,2)));
end
updated_para_PCIKF=ParY(:,end);
save updated_para_PCIKF.mat updated_para_PCIKF
total_iter_cost=sum(iter_num*2*(1+2*kl_term));
average_iter_cost=total_iter_cost/time_step;
%plot(RMSE);
%plot(Spread);








