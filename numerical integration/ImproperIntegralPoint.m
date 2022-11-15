function integralpoint = ImproperIntegralPoint(func,xlist)
while var(func(xlist))>1e-3
    xlist = xlist + (xlist(end)-xlist(1));
end
while length(xlist)>10 && var(func(xlist))<1e-3
    xlist = xlist(1:floor(length(xlist)/2));
end
integralpoint = xlist(end);
end