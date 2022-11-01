function [ integration ] = integralAdapt(func, left, right, length, prec)
%   自适应分段积分函数AdaptInt；
%   输入五个参数：被积函数（句柄）func，积分上下限right，left，要求精度prec，积分总长length；
%   输出一个参数：积分值integration；


if nargin == 4
    prec = 1e-6;
end

trapeInt = (right - left)*(func(left) + func(right))/2;
midInt = (right - left)*func((left + right)/2);
err = (real(trapeInt - midInt))/3;    % 由中点公式和梯形公式差估算误差
if (abs(err) < prec && (right - left) < length/5)   % 如果误差符合要求，则使用辛普森公式计算较精确结果
    integration = (right - left)*(func(left) + 4*func((left + right)/2) + func(right))/6;
else    % 否则，二分该段，分别进行自适应分段积分
    integration = integralAdapt(func, left, (left + right)/2, prec/2, length) + integralAdapt(func, (left + right)/2, right, prec/2, length);
end
end