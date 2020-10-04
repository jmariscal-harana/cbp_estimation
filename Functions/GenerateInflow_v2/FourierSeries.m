function F = FourierSeries(t, func, num)
% F = [a0  0
%      a1 b1
%      .. ..
%      an bn] num
% f(t) = a0/2 + sum{an·cos(n·w·t) + bn·sin(n·w·t)}
% f'(t) = sum{-an·n·w·sin(n·w·t) + bn·n·w·cos(n·w·t)}
    dt = t(2)-t(1);
    T = t(end);
    F = zeros(num + 1,2);
    for n = 0:num
        a0 = 0;
        an = 0;
        bn = 0;
        if n == 0
            for i = 2:length(func)
               a0 = a0 + (func(i)+func(i-1))/2*dt;
            end
            F(1,1) = 2/T*a0;
        else
            w = 2*pi/T;
            for i = 2:length(func)
                an = an + (func(i)*cos(n*w*t(i)) + func(i-1)*cos(n*w*t(i-1)))/2*dt; 
                bn = bn + (func(i)*sin(n*w*t(i)) + func(i-1)*sin(n*w*t(i-1)))/2*dt;
            end
            F(n + 1, 1) = 2/T*an;
            F(n + 1, 2) = 2/T*bn;
        end
    end
end