function integration = integralGlag(func,alpha,beta,n)

%------------- BEGIN CODE --------------

if nargin == 4
    [T,W] = gauss_laguerre(n);
else
    [T,W] = gauss_laguerre(5);
end

integration = 1/beta*exp(-alpha*beta)*dot(W,func(T/beta+alpha));

%------------- END OF CODE --------------
end

