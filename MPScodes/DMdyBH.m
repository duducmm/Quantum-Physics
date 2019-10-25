clear
clc

N=21;
Dmax=20;
Dwarm=10;
precision=10^(-12);
Epsilon=10^(-12);

Jx=1;
Jy=Jx;
Jz=0;
Un=1;
F=0.1;
h=2*Jx;
Wd=0;
Gamma=1;

dtwarm=0.05;
dt=0.01;

Ndtwarm=round(N/dtwarm)/2;
Ndt=round(N/dt)/2;
OsetZ=cell(N,N);
d=4;
dl=d^2;

a=zeros(d);
for I=2:d
    a(I-1,I)=(I-1)^(1/2);%anihilation operator
end

sx=(a+a')/2;
sy=1i*(a-a')/2;
sz=a'*a;
NLO=a'*sz*a;

Id=eye(dl);
id=eye(d);
for j=1:N,
    for jj=1:N,
        OsetZ{j,jj}=Id;
    end
end
Oseta=OsetZ;
for j=1:N
    OsetZ{j,j}=kron(id,sz);
    Oseta{j,j}=kron(id,a);
end
Fid=ones(1,Ndtwarm+Ndt);

mps=cell(1,N);
for j=1:N
    state=id(:,1);
    state=state*state';
    state=state/(trace(state));
    state=reshape(state,dl,1);
    mps{j}=reshape(state,[1,1,dl]);
end
% %mps=createrandommps(N,Dmin,dl);
mps=prepare(mps);

hz= Wd.*(rand(N,1)-0.5)*2+h;


% time evolution operator


H1=-1i*hz(1)*kron(kron(id,sz),Id)/2+1i*hz(1)*kron(kron(sz,id),Id)/2 ...
    -1i*F*kron(kron(id,sx),Id)/2+1i*F*kron(kron(sx,id),Id)/2 ...
    -1i*Un*kron(kron(id,NLO),Id)/2+1i*Un*kron(kron(NLO,id),Id)/2;
HN=-1i*hz(N)*kron(Id,kron(id,sz))/2+1i*hz(N)*kron(Id,kron(sz,id))/2 ...
    -1i*F*kron(Id,kron(id,sx))/2+1i*F*kron(Id,kron(sx,id))/2 ...
    -1i*Un*kron(Id,kron(id,NLO))/2+1i*Un*kron(Id,kron(NLO,id))/2;


J=sqrt(Gamma)*a;
L1=kron(kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)),Id)/2;
LN=kron(Id,kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)))/2;
Lj=L1+LN;

I=reshape(Id,[1,1,dl,dl]);
mpo_even=cell(1,N);
mpo_odd=cell(1,N);
for j=1:N
    mpo_even{j}=I;
    mpo_odd{j}=I;
end


for j=1:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2 ...
        +F*kron(kron(id,sx),Id)/2+F*kron(Id,kron(id,sx))/2 ...
        +Un*kron(kron(id,NLO),Id)/2+Un*kron(Id,kron(id,NLO))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2 ...
        +F*kron(kron(sx,id),Id)/2+F*kron(Id,kron(sx,id))/2 ...
        +Un*kron(kron(NLO,id),Id)/2+Un*kron(Id,kron(NLO,id))/2;
    LH=LH+1i*H+Lj;
    
    if j==1
        LH=LH+H1+L1;
    elseif j==(N-1)
        LH=LH+HN+LN;
    end
    w=expm(dtwarm*LH);
    
    [U,V] = BondMPOsvd(w,dl);
    
    mpo_odd{j}=U;
    mpo_odd{j+1}=V;
    
end
for j=2:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2 ...
        +F*kron(kron(id,sx),Id)/2+F*kron(Id,kron(id,sx))/2 ...
        +Un*kron(kron(id,NLO),Id)/2+Un*kron(Id,kron(id,NLO))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2 ...
        +F*kron(kron(sx,id),Id)/2+F*kron(Id,kron(sx,id))/2 ...
        +Un*kron(kron(NLO,id),Id)/2+Un*kron(Id,kron(NLO,id))/2;
    LH=LH+1i*H+Lj;
    
    if j==(N-1)
        LH=LH+HN+LN;
    end
    w=expm(dtwarm/2*LH);
    
    [U,V] = BondMPOsvd(w,dl);
    
    mpo_even{j}=U;
    mpo_even{j+1}=V;
end

% time evolution
TT=[linspace(1,Ndtwarm,Ndtwarm)*dtwarm Ndtwarm*dtwarm+linspace(1,Ndt,Ndt)*dt];
mpsSVD=mps;
for step=1:Ndtwarm
    clc
    step
    mpsi=mpsSVD;
    mpsSVD=MpoMps(mpo_even,mpsSVD,N);
    [mpsSVD]=mpsTruncSvd2(mpsSVD,Dwarm,Epsilon);
    mpsSVD=MpoMps(mpo_odd,mpsSVD,N);
    [mpsSVD]=mpsTruncSvd(mpsSVD,Dwarm,Epsilon);
    mpsSVD=MpoMps(mpo_even,mpsSVD,N);
    [mpsSVD]=mpsTruncSvd2(mpsSVD,Dwarm,Epsilon);
    mpsSVD=prepare(mpsSVD);
    
    Fid(step)=abs(braket(mpsi,mpsSVD))^2;
    
    [VZ,e,n]=expectationvalueTr(mpsSVD,OsetZ);
    
    Zsurf(:,step)=VZ.'/n;
end

mps=mpsSVD;

for j=1:N
    mpo_even{j}=I;
    mpo_odd{j}=I;
end


for j=1:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2 ...
        +F*kron(kron(id,sx),Id)/2+F*kron(Id,kron(id,sx))/2 ...
        +Un*kron(kron(id,NLO),Id)/2+Un*kron(Id,kron(id,NLO))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2 ...
        +F*kron(kron(sx,id),Id)/2+F*kron(Id,kron(sx,id))/2 ...
        +Un*kron(kron(NLO,id),Id)/2+Un*kron(Id,kron(NLO,id))/2;
    LH=LH+1i*H+Lj;
    
    if j==1
        LH=LH+H1+L1;
    elseif j==(N-1)
        LH=LH+HN+LN;
    end
    w=expm(dt*LH);
    
    [U,V] = BondMPOsvd(w,dl);
    
    mpo_odd{j}=U;
    mpo_odd{j+1}=V;
    
end
for j=2:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2 ...
        +F*kron(kron(id,sx),Id)/2+F*kron(Id,kron(id,sx))/2 ...
        +Un*kron(kron(id,NLO),Id)/2+Un*kron(Id,kron(id,NLO))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2 ...
        +F*kron(kron(sx,id),Id)/2+F*kron(Id,kron(sx,id))/2 ...
        +Un*kron(kron(NLO,id),Id)/2+Un*kron(Id,kron(NLO,id))/2;
    LH=LH+1i*H+Lj;
    
    if j==(N-1)
        LH=LH+HN+LN;
    end
    w=expm(dt/2*LH);
    
    [U,V] = BondMPOsvd(w,dl);
    
    mpo_even{j}=U;
    mpo_even{j+1}=V;
end

for step=1:Ndt
    
    clc
    step
    mpsi=mps;
    mps=MpoMps(mpo_even,mps,N);
    [mpsSVD]=mpsTruncSvd2(mps,Dmax,Epsilon);
    mpsSVD=prepare(mpsSVD);
    [mps]=VarCompressMPS(mps,mpsSVD,precision);
    
    mps=MpoMps(mpo_odd,mps,N);
    [mpsSVD]=mpsTruncSvd2(mps,Dmax,Epsilon);
    mpsSVD=prepare(mpsSVD);
    [mps]=VarCompressMPS(mps,mpsSVD,precision);
    
    mps=MpoMps(mpo_even,mps,N);
    [mpsSVD]=mpsTruncSvd2(mps,Dmax,Epsilon);
    mpsSVD=prepare(mpsSVD);
    [mps,K]=VarCompressMPS(mps,mpsSVD,precision);
    
    Fid(Ndtwarm+step)=abs(braket(mpsi,mps))^2;
    
    [VZ,e,n]=expectationvalueTr(mps,OsetZ);
    
    Zsurf(:,Ndtwarm+step)=VZ.'/n;
end

[VZ]=expectationvalueTr(mps,OsetZ);
[Va,e,n]=expectationvalueTr(mps,Oseta);

for nn=1:N
    for j=1:N,
        osetg2{1,j}=Id;
    end
    osetg2{1,ceil(N/2)}=kron(id,sz);
    osetg2{1,nn}=kron(id,sz);
    if nn==ceil(N/2)
        osetg2{1,nn}=kron(id,NLO);
    end
    [~,e,Tr]=expectationvalueTr(mps,osetg2);
    Spacialg2(nn)=e/Tr;
end
Spacialg2=Spacialg2'./(VZ.'/n*VZ(ceil(N/2)).'/n);
for nn=1:N
    for j=1:N,
        osetg1{1,j}=Id;
    end
    osetg1{1,ceil(N/2)}=kron(id,a');
    osetg1{1,nn}=kron(id,a);
    if nn==ceil(N/2)
        osetg1{1,nn}=kron(id,sz);
    end
    [~,e,Tr]=expectationvalueTr(mps,osetg1);
    Spacialg1(nn)=e/Tr;
end
Spacialg1=abs(Spacialg1)'./(abs(Va.'/n*Va(ceil(N/2))'/n));

figure(2),
%pcolor(TT,[1:1:N],(real(Zsurf))), shading interp
 plot(TT,log10(1-Fid))

hold on
subplot(131),plot([1:1:N],real(VZ.'/n),'o-')
hold on
subplot(132),plot([1:1:N],real(Spacialg1),'o-')
hold on
subplot(133),plot([1:1:N],real(Spacialg2),'o-')

hold on
plot([0:1:floor(N/2)],real(Spacialg2(ceil(N/2):end)),'o-')

for rea=1:Nrea
Vh(rea)=2*Jx-(rea-1)/(Nrea-1)*4*Jx;
end
pcolor([1:1:N],-Vh,real(Mg2)'), shading interp
pcolor([1:1:N],-Vh,real(Mg1)'), shading interp
pcolor([1:1:N],-Vh,real(Mvz)'), shading interp

