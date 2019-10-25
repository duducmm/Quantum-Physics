function [Q,R]=qr2(T)
[m,n]=size(T);
Delt=10^(-13);
    try
    [Q,R]=qr(T,0);
    catch
        [Q,R]=qr(T+Delt*randn(m,n),0);
    end