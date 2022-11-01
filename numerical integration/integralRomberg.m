function integration = integralRomberg(func, left, right)

prec = 1e-5;

N = 2;
T = zeros(N,N);
w = right-left;

T(1,1) = w/2*(func(left)+func(right));
i = 2;
j = 1;
while 1
    while i>=j
        if j == 1
            w0 = w/2^(i-1);
            xlst = left:w0:right;
            T(i,j) = w0/2*sum(func(xlst(1:end-1))+func(xlst(2:end)));
        else
            T(i,j) = (4*T(i,j-1)-T(i-1,j-1))/3;
        end
        if j == size(T,1)
            Ttemp = zeros(j+1,j+1);
            Ttemp(1:j,1:j) = T;
            T = Ttemp;
        end
        j = j+1;
    end
    a = T(end-1,end-1);
    if abs(T(i,i)-T(i-1,i-1))<prec && j>4
        break;
    end
    j = 1;
    i = i + 1;
end
integration = T(end-1,end-1);
end