function integration = integralGLOfPoly(objfun,fullfun,xList,n)

if nargin == 4
    [T,W] = gauss_lege(n);
else
    [T,W] = gauss_lege(5);
end

% exp_func = @(x) exp(1i.*abs(POSITION.Z).*kz(x));
yList = objfun(xList);
[~,~,~,~,C] = cubicSplineInterpolation(xList, yList);
integration = 0;
for ii = 1: length(xList)-1
    func = @(x)C(ii,1)+C(ii,2).*x+C(ii,3).*x.^2+C(ii,4).*x.^3;
    func = @(x)func(x).*fullfun(x)./objfun(x);
    leftEndpoint = xList(ii);
    rightEndpoint = xList(ii+1);
    xlst = (leftEndpoint+rightEndpoint)/2+(rightEndpoint-leftEndpoint)/2*T;
    integration = integration + (rightEndpoint-leftEndpoint)/2*dot(W,func(xlst));
end