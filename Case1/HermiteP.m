function y = HermiteP(x)
% 1D Hermite Polynomial
y = zeros(3,1);
y(1) = 1;
y(2) = x;
y(3) = (x^2-1)/sqrt(2);
end