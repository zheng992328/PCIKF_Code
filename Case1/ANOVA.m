function F1 = ANOVA(CP,C0,C1,NR)
% map the collocation point to model parameter through 1st order Anova
% approximation

% CP: collocation point
% C0: 0th order coefficient
% C1: 1st order coefficient
F1 = C0;
for i=1:NR
    t = HermiteP(CP(i));
    F1 = F1+C1(:,:,i)*t(2:3);
end