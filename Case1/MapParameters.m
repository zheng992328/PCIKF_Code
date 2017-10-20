function p = MapParameters(si,NR, C_P0, C_P1)
% mapping from random space to parameter space, using PCE representation of
% the parameters

p = C_P0;
HP = zeros(3,NR); 
for i = 1:NR
    HP(:,i) = HermiteP(si(i));
    p = p + C_P1(:,:,i)*HP(2:3,i);
end
% for i = 1:NR
%     for j = 1:NR
% % for j = 1:NR-1
% %     for i = j+1:NR
%         for u =1:2
%             for v =1:2
%                 p = p+C_P2(:,v,u,j,i)*HP(1+v,j)*HP(1+u,i);
%             end
%         end
%     end
% end
