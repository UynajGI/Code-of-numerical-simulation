function integration = integralGL(xList,func,n)

%------------- BEGIN CODE --------------

if nargin == 3
    [T,W] = gauss_lege(n);
else
    [T,W] = gauss_lege(5);
end

X_LENGTH = length(xList);
integration = 0;

for i = 1:X_LENGTH-1
    leftEndpoint = xList(i);
    rightEndpoint = xList(i+1);
    xlst = (leftEndpoint+rightEndpoint)/2+(rightEndpoint-leftEndpoint)/2*T;
    integration = integration + (rightEndpoint-leftEndpoint)/2*dot(W,func(xlst));
end

%------------- END OF CODE --------------
end

