function y = ogs_call(x,i,t)
% save (strcat('.\gas_',num2str(i),'\k_H.mat'),'x');
cd (strcat('./gas_',num2str(i)));
save('x.mat','x');
system('python ./write_mmp_multi_para.py');
system(strcat('python ./time_modify.py ',32,num2str(t)));
system('./ogs  first_case >>a 2>err');
% pause(12)
system('python ./read_values.py ');
data=importdata(strcat('data_all_time_',num2str(i),'.txt'));
y=data; 
cd ('..')
end
 
