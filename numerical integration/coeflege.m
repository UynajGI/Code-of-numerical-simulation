function [a,b] = coeflege(n)
a = zeros(n,1);
b = a;
b(1) = 2;
k = 2:n;
b(k) = 1./(4-1./(k-1).^2);
end