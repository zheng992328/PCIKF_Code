function [C_O0,C_O1] = Get_ANOVA_adaptive(cur_path,kl_term,Node_Num,C_P0,C_P1,tt,Obs_Node,Obs_Num)
% global Npair C_O2

CP0 = zeros(kl_term,1);   % the 0th order collocation point is the anchor point

% cd(strcat(cur_path,'\start\start (1)'))

p = MapParameters(CP0,kl_term,C_P0,C_P1);  %%将CP0映射成可以代入原始模型的参数

% C_O0 = Hydrus(p,Node_Num,tt,Obs_Node);
C_O0=ogs_call(p,0,tt);

% 1st order ANOVA
P = [-sqrt(3); 0; sqrt(3)]; % 1D collocation points 
PC1 = zeros(3); % 1D polynomial chaos evaluated at the 3 collocation points
for i =1:3
    PC1(:,i) = HermiteP(P(i));
end

C_O1 = zeros(Obs_Num,2,kl_term);
C_O0_a = zeros(Obs_Num,1); % amendment to C_O0
% A matrix
A1 = zeros(3);
for j=1:3
    A1(j,:) = PC1(:,j);
end

EnergyO = zeros(Obs_Num,2,kl_term);  % normalized energy caused by each active dimension linear part and non linear part

parfor i = 1:kl_term
    fprintf('1_st_%d_1\n',i);
%     cd(strcat(cur_path,'\gas_ (',num2str(i),')'))
    b = zeros(3,Obs_Num);        % RHS vector b
    CP = zeros(kl_term,1); CP(i) = P(1);    
    p = MapParameters(CP,kl_term,C_P0,C_P1);
%     b(1,:) = (Hydrus(p,Node_Num,tt,Obs_Node)-C_O0); 
    b(1,:)=ogs_call(p,i,tt)-C_O0;

    fprintf('1_st_%d_2\n',i);
    CP = zeros(kl_term,1); CP(i) = P(3);
    p = MapParameters(CP,kl_term,C_P0,C_P1);
%     b(3,:) = (Hydrus(p,Node_Num,tt,Obs_Node)-C_O0);
    b(3,:)=ogs_call(p,i,tt)-C_O0;
    
    k=A1\b;
    C_O1(:,:,i)=k(2:3,:)';
    C_O0_a = C_O0_a + k(1,:)';
    EnergyO(:,:,i) = C_O1(:,:,i).*C_O1(:,:,i);
end

C_O0 = C_O0 + C_O0_a;
% % 2nd order ANOVA
% C_O2 = zeros(Obs_Num,2,2,kl_term,kl_term);
% C_O0_t = C_O0;        % temp C_O0
% C_O1_t = C_O1;  % temp C_O1
% 
% % A matrix
% A2 = zeros (9);
% n=0;
% for u =1:3
%     for v =1:3
%         t = PC1(:,v)*PC1(:,u)';
%         n=n+1;
%         A2(n,:)=reshape(t,1,9);
%     end
% end
% 
% % adaptive criteria for 2nd order ANOVA
% Pairs = zeros(2,0);
% Npair = 0;
% % Energy_1st_t=sum(sum(EnergyO,2),3);  % total 1st order ANOVA energy of each output variables 
% Energy_1st=reshape(sum(EnergyO,2),Obs_Num,kl_term);
% N_Total_pair = kl_term*(kl_term-1)/2;
% Thresh = 0.05;    % 2nd order criterion, !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
% 
% for iO=1:Obs_Num
% %     break %  add this line to have 1st order ANOVA 
%     Pair_all = zeros(N_Total_pair,4);
%     n=0;
%     
%     for u=1:kl_term-1
%         for v=u+1:kl_term
%             n=n+1;
%             Pair_all(n,:)=[Energy_1st(iO,u)*Energy_1st(iO,v);0;u;v];
%         end
%     end
%     Pair_all = sortrows(Pair_all,-1);
%     
%     MaxPair = 50;  % control it by yourself
%     
%     for iPair=1:MaxPair
%         % compute the 2nd order ANOVA and estimate the error of estimation of Var(Output(j))
%         [Err_2nd,C_O0_t, C_O1_t, Pairs, Pair_all] = compute_2nd_ANOVA(cur_path,Pair_all, iPair, C_O0_t, C_O1_t, A2, Pairs,iO,C_O0,C_O1,kl_term,Node_Num,C_P0,C_P1,C_P2,tt,Obs_Node,Obs_Num); 
%         Total_out_E = sum(sum(C_O1_t(iO,:,:).^2))+sum(sum(sum(sum(C_O2(iO,:,:,:,:).^2))));       
% %         str = [num2str(iO) '  ' num2str(MaxPair) '  ' num2str(iPair) '\n'];
% %         fprintf(str);     
%         if Err_2nd < Thresh*(Total_out_E+1)       
% %         if Err_2nd < Thresh*Total_out_E
%             break
%         end
%     end
%     
% end
% 
% C_O0 = C_O0_t;    % update 0th and 1st order ANOVA terms AFTER all 2nd ANOVA terms calculated
% C_O1 = C_O1_t; 
% 
% 
