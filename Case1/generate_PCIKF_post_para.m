
function [para_sample,para_mean]=generate_PCIKF_post_para(t,N)
    kl_term=15;
    C_P=importdata(['./updated_para_sample_PCIKF_',num2str(t),'.mat']); %导入C_P
    C_P0=importdata(['./updated_para_mean_PCIKF_',num2str(t),'.mat']);   %导入C_P0
    C_P_1=C_P(:,2:end);
    C_P1=reshape(C_P_1,kl_term,2,kl_term);
    para_sample=zeros(kl_term,N);
    
    for i=1:N
        si=randn(kl_term,1);
        para_sample(:,i)=MapParameters(si,kl_term,C_P0,C_P1);
        
    end
	para_mean=C_P0;
    
    
    