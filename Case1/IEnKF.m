function [rmse,spread,average_iter_cost]=IEnKF(index,N)   
	% N=80;
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
	para_mean=reshape(para_mean,length(para_mean),1);
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
	load para_true.mat
	obs=ogs_call(para_true,999,time_step);
	obs_mat=zeros(time_step,Obs_Num,N);
	for i=1:time_step
		obs_mat(i,:,:)=repmat(obs((i-1)*Obs_Num+1:i*Obs_Num),1,N)+sqrt(varR)*randn(Obs_Num,N);
	end

	%% 初始参数
	para_initial=repmat(para_mean,1,N)+Fi*randn(kl_term,N);
	para_initial_ave=mean(para_initial,2);
	m2=para_initial;
	m1=para_initial;
	%% 反演部分
	Num=zeros(1,time_step);
	ParY=zeros(kl_term,time_step);%放更新之后的参数场
	RMSE=zeros(1,time_step);
	Spread=zeros(1,time_step);
	total_Imax=zeros(1,time_step);
	iter_num=zeros(1,time_step);

	for t=1:time_step
		fprintf('this is the %d time step\n',t);
		mpr=m2;
		para_error=mpr-repmat(mean(mpr,2),1,N);
		Cm=para_error*para_error'/(N-1);
		obs_mat_t=squeeze(obs_mat(t,:,:));
		obs_mat_t_error=obs_mat_t-repmat(mean(obs_mat_t,2),1,N);
		Cd=obs_mat_t_error*obs_mat_t_error'/(N-1);
		m1=mpr;
		Imax=1;
		while 1
			fprintf('this is the %d time %d iteration\n',t,Imax)
			m1_error=m1-repmat(mean(m1,2),1,N);
			obs_prediction=zeros(Obs_Num,N);
			parfor i=1:N
				obs_prediction(:,i)=ogs_call(m1(:,i),i,t);
			end
			obs_prediction_error=obs_prediction-repmat(mean(obs_prediction,2),1,N);
			
			wn_M=pinv(m1_error);
			G=obs_prediction_error*wn_M;
			m2=beta*mpr+(1-beta)*m1-beta*Cm*G'*pinv(Cd+G*Cm*G')*(obs_prediction-obs_mat_t-G*(m1-mpr));
			S1=trace(((obs_prediction-obs_mat_t)')*pinv(Cd)*(obs_prediction-obs_mat_t));
			
			parfor i=1:N
				obs_prediction(:,i)=ogs_call(m2(:,i),i,t);
			end
			S2=trace(((obs_prediction-obs_mat_t)')*pinv(Cd)*(obs_prediction-obs_mat_t));
			fprintf('S1: %d\n',S1);
			fprintf('S2: %d\n',S2);
			fprintf('max(m2-m1): %d\n',max(max(abs(m2-m1))));
			fprintf('S2-S1: %d\n',S2-S1);
			fprintf('eps2*S1: %d\n',eps2*S1);
			fprintf('\n');
			 if S2<S1
				if max(max(abs(m2-m1)))<eps1 | S1-S2<eps2*S1 | Imax>2
					break
				end
				m1=m2;
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
		ParY(:,t)=mean(m2,2);
		RMSE(t)=sqrt(mean((ParY(:,t)-para_true).^2));
		Spread(t)=sqrt(mean(var(m2,0,2)));
	end
	updated_para_EnRML=ParY(:,end);
	save updated_para_EnRML.mat updated_para_EnRML
	save(['./updated_para_EnRML_',num2str(index),'_',num2str(N),'.mat'],'updated_para_EnRML')
	total_iter_cost=sum(iter_num*2*N);
	average_iter_cost=total_iter_cost/time_step;
	rmse=RMSE(end);
	spread=Spread(end);
	%plot(RMSE);
	%plot(Spread);


