total_index=10;
kl_term=15;
updated_para_IEnKF_80_mat=zeros(3,total_index);

for i=1:total_index
    fprintf('this is %d calculation\n',i);
    [rmse,spread,average_iter_cost]=IEnKF(i,80);
    updated_para_IEnKF_80_mat(1,i)=rmse;
    updated_para_IEnKF_80_mat(2,i)=spread;
    updated_para_IEnKF_80_mat(3,i)=average_iter_cost;
end

updated_IEnKF_80=mean(updated_para_IEnKF_80_mat,2)


updated_para_IEnKF_100_mat=zeros(3,total_index);
for i=1:total_index
    fprintf('this is %d calculation\n',i);
    [rmse,spread,average_iter_cost]=IEnKF(i,100);
    updated_para_IEnKF_100_mat(1,i)=rmse;
    updated_para_IEnKF_100_mat(2,i)=spread;
    updated_para_IEnKF_100_mat(3,i)=average_iter_cost;
end
updated_IEnKF_100=mean(updated_para_IEnKF_100_mat,2)


