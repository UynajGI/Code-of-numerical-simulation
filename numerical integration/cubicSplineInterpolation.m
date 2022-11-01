function [D,h,d,M,C] = cubicSplineInterpolation(xList, yList)


%------------- BEGIN CODE --------------

X_LENGTH = length(xList);
xList = reshape(xList,X_LENGTH,1);
yList = reshape(yList,X_LENGTH,1);
D=zeros(X_LENGTH,X_LENGTH);

h = diff(xList);
u = [0;h(1:end-1)./(h(1:end-1)+h(2:end));0];
n = 1-u;

D(1:X_LENGTH+1:end) = 2;
D(2:X_LENGTH+1:end) = u(2:end);
D(X_LENGTH+1:X_LENGTH+1:end) = n(1:end-1);


j = 2:X_LENGTH-1;
d = 6*(h(j-1).*(yList(j+1)-yList(j)))./(h(j).*(h(j)+h(j-1)))-6*(h(j).*(yList(j)-yList(j-1)))./(h(j-1).*(h(j)+h(j-1)));
d1 = d(1);
dN1 = d(end);
d = [d1;d;dN1];

M = D\d;

C = zeros(X_LENGTH-1,4);
j = 1:X_LENGTH-1;
C(:,1) = h(j)./6.*(xList(j).*M(j+1)-xList(j+1).*M(j))+1./(6.*h(j)).*(xList(j+1).^3.*M(j)-xList(j).^3.*M(j+1))-1./h(j).*(xList(j).*yList(j+1)-xList(j+1).*yList(j));
C(:,2) = h(j)./6.*(M(j)-M(j+1))+1./(2.*h(j)).*(-xList(j+1).^2.*M(j)+xList(j).^2.*M(j+1))+1./h(j).*(yList(j+1)-yList(j));
C(:,3) = 1./(2.*h(j)).*(xList(j+1).*M(j)-xList(j).*M(j+1));
C(:,4) = 1./(6.*h(j)).*(M(j+1)-M(j));

%------------- END OF CODE --------------
end