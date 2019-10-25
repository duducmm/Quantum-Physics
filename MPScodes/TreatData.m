clear
clc
VN=[20 25 30 35];
VWd=[0 0.5 1 1.5 2];
Gamma=1;
Jz=1;
mu=1;
Wd=0;
Dmax=50;

Vcurrent=[];
VErrCurrent=[];
gap=[];
options = optimset('TolFun',1e-10,'TolFun',1e-10);

for i=1:4

    N=VN(i);
fname = sprintf('DMdyJz%dW%dG%dmu%dDmax%dN%d.mat',Jz,Wd,Gamma,mu,Dmax,N);
    load(fname)
    
Vcurrent=[Vcurrent current];
%VErrCurrent=[VErrCurrent ErrCurrent];
VErrCurrent=[VErrCurrent Err];
hold on
%errorbar([0:1/(N-1):1],real(VZ.'/n),imag(VZ.'/n),'o-')
%errorbar(real((Vcurr(1:N-1)-Vcurr2(1:N-1))/n),imag((Vcurr(1:N-1)-Vcurr2(1:N-1))/n),'o-')
 plot(TT,(1-Fid))
% 
% [alpha ,Dist] = fminsearch(@(alpha)ffiitt((1-Fid(Ndtwarm+1:end)),TT(Ndtwarm+1:end),alpha),[-rand rand]...
%     ,options);
% plot(TT,alpha(2)*TT.^alpha(1))
% gap=[gap -alpha(1)];

end

plot(VN,gap)

hold on
errorbar(VN,Vcurrent,VErrCurrent,'o')

[alpha ,Dist] = fminsearch(@(alpha)ffiitt(Vcurrent,VN,alpha),[-rand rand]...
    ,options);

hold on
plot(VN,alpha(2)*VN.^alpha(1))


[alpha ,Dist] = fminsearch(@(alpha)ffiiexp(Vcurrent,VN,alpha),[-rand rand]...
    ,options);
hold on
plot(VN,abs(alpha(1)).^(VN*alpha(2)))

