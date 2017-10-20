total_index=10;

%% 计算参考的均值和标准差
load updated_para_sample_EnRML_1_500.mat
mu_b=mean(m2,2);
sigma_b=std(m2,0,2);

%% 计算PCIKF对应的H2
load updated_para_PCIKF.mat
load updated_para_sample_PCIKF
mu_PCIKF=updated_para_PCIKF;
sigma_PCIKF=(sum(C_P(:,2:end).^2,2)).^0.5;
h=zeros(size(mu_PCIKF));
for i=1:length(mu_PCIKF)
	expon=-((mu_PCIKF(i)-mu_b(i))^2)/(4*(sigma_PCIKF(i)^2+sigma_b(i)^2));
	dishu=sqrt(2*sigma_PCIKF(i)*sigma_b(i)/(sigma_PCIKF(i)^2+sigma_b(i)^2));
	h(i)=1-dishu*exp(expon);
end
H_PCIKF=mean(h)

%% 计算样本为30时候的礖2
H_30_=zeros(total_index,1);
for k=1:total_index
	load( ['updated_para_sample_EnRML_',num2str(k),'_',num2str(30),'.mat'])
	mu_30=mean(m2,2);
	sigma_30=std(m2,0,2);
	h=zeros(size(mu_30));
	for i=1:length(mu_30)
		expon=-((mu_30(i)-mu_b(i))^2)/(4*(sigma_30(i)^2+sigma_b(i)^2));
		dishu=sqrt(2*sigma_30(i)*sigma_b(i)/(sigma_30(i)^2+sigma_b(i)^2));
		h(i)=1-dishu*exp(expon);
	end
	H_30=mean(h);
	H_30_(k)=H_30;
end
H_30_average=mean(H_30_)


%% 计算样本为50时候的礖2
H_50_=zeros(total_index,1);
for k=1:total_index
	load( ['updated_para_sample_EnRML_',num2str(k),'_',num2str(50),'.mat'])
	mu_50=mean(m2,2);
	sigma_50=std(m2,0,2);
	h=zeros(size(mu_50));
	for i=1:length(mu_50)
		expon=-((mu_50(i)-mu_b(i))^2)/(4*(sigma_50(i)^2+sigma_b(i)^2));
		dishu=sqrt(2*sigma_50(i)*sigma_b(i)/(sigma_50(i)^2+sigma_b(i)^2));
		h(i)=1-dishu*exp(expon);
	end
	H_50=mean(h);
	H_50_(k)=H_50;
end
H_50_average=mean(H_50_)
