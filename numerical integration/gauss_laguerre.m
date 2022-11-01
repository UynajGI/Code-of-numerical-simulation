function [x,w] = gauss_laguerre(n)
% Nodes and weights for Gauss-Laguerre quadrature of
% arbitrary order
% Input: n =numbers of nodes of quadrature rule
% Output: x= vector of nodes, w = vector of weights
J = diag(1:2:2*n-1)- diag(1:n-1,1) - diag(1:n-1,-1);
% Jacobi matrix
[V,D] = eig(J); 
[x,ix] = sort(diag(D));
% nodes are eigenvalues, which are on diagonal of D
w= 1*V(1,ix)'.^2;
% V(1,ix)’ is column vector of first row of sorted V
end