function integration = integralcubicSpline(xList, yList)

%------------- BEGIN CODE --------------

[~,~,~,~,C] = cubicSplineInterpolation(xList, yList);

X_LENGTH=length(xList); 
xList = reshape(xList,X_LENGTH,1);
integration = trace(C*(1./(1:4).*(xList(2:end).^(1:4)-xList(1:end-1).^(1:4)))');

%------------- END OF CODE --------------
end