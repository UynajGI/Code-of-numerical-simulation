function resultsList = fitcubicSpline(xList, yList, fitList)


%------------- BEGIN CODE --------------

[~,~,~,~,C] = cubicSplineInterpolation(xList, yList);

fitList = sort(fitList);
FIT_LENGTH = length(fitList);
resultsList = zeros(1,FIT_LENGTH);

istrat = 1;

for ifit = 1:FIT_LENGTH
    for i = istrat:X_LENGTH-1
        if xList(i)<=fitList(ifit)&&xList(i+1)>=fitList(ifit)
            resultsList(ifit) = dot(C(i,:),(fitList(ifit).^(0:3)));
            istrat = i;
            break
        end
    end
end

%------------- END OF CODE --------------
end