function [Err_2nd, C_O0_t, C_O1_t, Pairs, Pair_all] = compute_2nd_ANOVA(cur_path,Pair_all, iPair,C_O0_t, C_O1_t, A2, Pairs, iO, C_O0, C_O1, kl_term,Node_Num, C_P0, C_P1,C_P2,tt,Obs_Node,Obs_Num)
                                                                       
% compute the 2nd order ANOVA and estimate the error of estimation of Var(Output(j))
global Npair C_O2
    newPair = 1;
    for i=1:Npair
        if max(abs(Pair_all(iPair,3:4)-Pairs(:,i)'))==0
            newPair = 0;
            break
        end
    end
    
    if newPair
        Npair=Npair+1;
    
        P = [-sqrt(3); 0; sqrt(3)]; % 1D collocation points 

        % right hand side vector b
        
        CP = zeros(kl_term,1); CP(Pair_all(iPair,3))=P(1); CP(Pair_all(iPair,4))=P(1);
        CP1(:,1) = CP;
        CP = zeros(kl_term,1); CP(Pair_all(iPair,3))=P(3); CP(Pair_all(iPair,4))=P(1);
        CP1(:,2) = CP;
        CP = zeros(kl_term,1); CP(Pair_all(iPair,3))=P(1); CP(Pair_all(iPair,4))=P(3);
        CP1(:,3) = CP;
        CP = zeros(kl_term,1); CP(Pair_all(iPair,3))=P(3); CP(Pair_all(iPair,4))=P(3);
        CP1(:,4) = CP;
        save 'CP1.mat' CP1;
        save 'C_O0.mat' C_O0
        save 'C_O1.mat' C_O1
        
        D=zeros(4,Obs_Num);
        for i = 1:4
            fprintf('2_nd_%d_%d\n',iO,i);
%             cd(strcat(cur_path,'\start\start (',num2str(i),')'))
            p= MapParameters(CP1(:,i),kl_term,C_P0,C_P1,C_P2);
%             D(i,:) = (Hydrus(p,Node_Num,tt,Obs_Node)-ANOVA(CP1(:,i),C_O0,C_O1,kl_term))';
            zz=ogs_call(p,i,tt);
            D(i,:)=zz-ANOVA(CP1(:,i),C_O0,C_O1,kl_term);             
        end
        b = zeros(9,Obs_Num);
        b(1,:) = D(1,:);
        b(3,:) = D(2,:);
        b(7,:) = D(3,:);
        b(9,:) = D(4,:);
         
        k=A2\b;
        for s =1:Obs_Num
            t=reshape(k(:,s),3,3);
            C_O0_t(s) = C_O0_t(s)+t(1,1);
            C_O1_t(s,:,Pair_all(iPair,4)) = C_O1_t(s,:,Pair_all(iPair,4))+t(1,2:3);
            C_O1_t(s,:,Pair_all(iPair,3)) = C_O1_t(s,:,Pair_all(iPair,3))+t(2:3,1)';
            C_O2(s,:,:,Pair_all(iPair,3),Pair_all(iPair,4)) = t(2:3,2:3);
            if s == iO
                Pair_all(iPair,2)=sum(sum(t(2:3,2:3).^2));
            end
        end
        
        % Ratio = mean(Pair_all(1:iPair,2)./(Pair_all(1:iPair,1)+(1e-10)));
        Ratio= sum(Pair_all(1:iPair,2))/(sum(Pair_all(1:iPair,1))+(1e-10));
        Err_2nd = Ratio*sum(Pair_all(iPair+1:end,1));
        Pairs = [Pairs,Pair_all(iPair,3:4)'];

    else
        % if not new pair
        Pair_all(iPair,2)=sum(reshape(C_O2(iO,:,:,Pair_all(iPair,3),Pair_all(iPair,4)).^2,4,1));
        % Ratio = mean(Pair_all(1:iPair,2)./(Pair_all(1:iPair,1)+(1e-10)));
        Ratio= sum(Pair_all(1:iPair,2))/(sum(Pair_all(1:iPair,1))+(1e-10));
        Err_2nd = Ratio*sum(Pair_all(iPair+1:end,1));
    end