function [x,w] = gauss_lege(n)
[a,b] = coeflege(n);
JacM = diag(a) + diag(sqrt(b(2:n)),1) + diag(sqrt(b(2:n)),-1);
[w,x] = eig(JacM);
x = diag(x);
scal = 2;
w = w(1,:)'.^2*scal;
[x,ind] = sort(x);
w = w(ind);
end