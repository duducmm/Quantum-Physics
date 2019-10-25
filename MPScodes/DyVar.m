clear
clc

N=20;
Dwarm=20;
Dmax=50;
precision=10^(-10);
Epsilon=10^(-10);

Jx=1;
Jy=1;
Jz=1;
h=0;
Wd=0;
NsweepsMax=8;
Gamma=2;
GAMMA=zeros(N,1);
GAMMA(1)=Gamma;
GAMMA(N)=Gamma;


dtwarm=0.1;
dt=0.05;

Ndtwarm=round(2*N/dtwarm);
Ndt=round(2*N/dt);
Fid=ones(1,Ndtwarm+Ndt);
% magnetization in z-direction
OsetZ=cell(N,N);
OsetCurrent=cell(N-1,N);
d=2;
dl=d^2;
sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1];
s=[0 0 ;1 0];
Id=eye(dl);
id=eye(d);
for j=1:N,
    for jj=1:N,
        OsetZ{j,jj}=Id;
        OsetCurrent{j,jj}=Id;
    end
end
OsetCurrent2=OsetCurrent;
for j=1:N-1
    OsetZ{j,j}=kron(id,sz);
    OsetCurrent{j,j}=kron(id,sx);
    OsetCurrent{j,j+1}=kron(id,sy);
    OsetCurrent2{j,j}=kron(id,sy);
    OsetCurrent2{j,j+1}=kron(id,sx);
end
OsetZ{N,N}=kron(id,sz);

% starting state
mps=cell(1,N);
for j=1:N
    %state=[1; 0];
    %state=state*state';
    state=eye(d);
    state=state/(trace(state));
    state=reshape(state,dl,1);
    mps{j}=reshape(state,[1,1,dl]);
end
% %mps=createrandommps(N,Dmin,dl);
mps=prepare(mps);

hz= Wd.*(rand(N,1)-0.5)*2+h;


% time evolution operator


H1=-1i*hz(1)*kron(kron(id,sz),Id)/2+1i*hz(1)*kron(kron(sz,id),Id)/2;
HN=-1i*hz(N)*kron(Id,kron(id,sz))/2+1i*hz(N)*kron(Id,kron(sz,id))/2;

J=sqrt(Gamma)*s';
L1=kron(kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)),Id);
J=sqrt(Gamma)*s;
LN=kron(Id,kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)));

%Llocal=Llocal+kron(L,Id)/2+kron(Id,L)/2;




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
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2;
    LH=LH+1i*H;
    
    w=expm(dtwarm*LH);
    
    if j==1
        w=expm(dtwarm*(LH+H1+L1));
    elseif j==(N-1)
        w=expm(dtwarm*(LH+HN+LN));
    end
    
    [U,V] = BondMPOsvd(w,dl);
    
    mpo_odd{j}=U;
    mpo_odd{j+1}=V;
    
end
for j=2:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2;
    LH=LH+1i*H;
    
    w=expm(dtwarm/2*LH);
    
    if j==(N-1)
        w=expm(dtwarm/2*(LH+HN+LN));
        
    end
    
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


for j=1:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2;
    LH=LH+1i*H;
    
    w=expm(dt*LH);
    
    if j==1
        w=expm(dt*(LH+H1+L1));
    elseif j==(N-1)
        w=expm(dt*(LH+HN+LN));
    end
    
    [U,V] = BondMPOsvd(w,dl);
    
    mpo_odd{j}=U;
    mpo_odd{j+1}=V;
    
end
for j=2:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2;
    LH=LH+1i*H;
    
    w=expm(dt/2*LH);
    
    if j==(N-1)
        w=expm(dt/2*(LH+HN+LN));
        
    end
    
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
[Vcurr,Tcur,]=expectationvalueTr(mps,OsetCurrent);
[Vcurr2,Tcur2,n]=expectationvalueTr(mps,OsetCurrent2);
current=real((Tcur-Tcur2)/(N-1)/n);
Err=imag((Tcur-Tcur2)/(N-1)/n);


M=N+3*(N-2)+3+3*(N-2)+3;
Lset=cell(M,N);
J=[0 0; 1 0]';

for m=1:M, 
    for j=1:N, 
        Lset{m,j}=Id;               
    end
end


for j=1:(N-1)
Lset{N+3*(j-1)+1,j}=-1i*sqrt(Jx)*kron(id,sx);
Lset{N+3*(j-1)+1,j+1}=sqrt(Jx)*kron(id,sx);
Lset{N+3*(j-1)+2,j}=-1i*sqrt(Jy)*kron(id,sy);
Lset{N+3*(j-1)+2,j+1}=sqrt(Jy)*kron(id,sy);
Lset{N+3*(j-1)+3,j}=-1i*sqrt(Jz)*kron(id,sz);
Lset{N+3*(j-1)+3,j+1}=sqrt(Jz)*kron(id,sz);

Lset{N+3*(N-2)+3+3*(j-1)+1,j}=1i*sqrt(Jx)*kron(sx.',id);
Lset{N+3*(N-2)+3+3*(j-1)+1,j+1}=sqrt(Jx)*kron(sx.',id);
Lset{N+3*(N-2)+3+3*(j-1)+2,j}=1i*sqrt(Jy)*kron(sy.',id);
Lset{N+3*(N-2)+3+3*(j-1)+2,j+1}=sqrt(Jy)*kron(sy.',id);
Lset{N+3*(N-2)+3+3*(j-1)+3,j}=1i*sqrt(Jz)*kron(sz.',id);
Lset{N+3*(N-2)+3+3*(j-1)+3,j+1}=sqrt(Jz)*kron(sz.',id);

% Lset{N+1*(j-1)+1,j}=-1i*sqrt(Jx)*kron(id,sx);
% Lset{N+1*(j-1)+1,j+1}=sqrt(Jx)*kron(id,sx);
% Lset{N+1*(N-2)+1+1*(j-1)+1,j}=1i*sqrt(Jx)*kron(sx.',id);
% Lset{N+1*(N-2)+1+1*(j-1)+1,j+1}=sqrt(Jx)*kron(sx.',id);

Lset{j,j}=GAMMA(j)*(kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)))...
    -1i*hz(j)*kron(id,sz)+1i*hz(j)*kron(sz.',id);

%Lset{N+1,j}=Vec0*reshape(id,1,2^2);

end
j=j+1;
J=[0 0; 1 0];
Lset{j,j}=GAMMA(j)*(kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)))...
     -1i*hz(j)*kron(id,sz)+1i*hz(j)*kron(sz.',id);

% for j=1:(N-2)
% Lset{N+3*(N-2)+3+3*(N-2)+3+(j-1)+1,j}=-1i*sqrt(Jx2)*kron(id,sx);
% Lset{N+3*(N-2)+3+3*(N-2)+3+(j-1)+1,j+2}=sqrt(Jx2)*kron(id,sx);
% Lset{N+3*(N-2)+3+3*(N-2)+3+N-2+(j-1)+1,j}=1i*sqrt(Jx2)*kron(sx.',id);
% Lset{N+3*(N-2)+3+3*(N-2)+3+N-2+(j-1)+1,j+2}=sqrt(Jx2)*kron(sx.',id);
% end
 

% Asymptotic State

[E0,mps0]=minimizeL(Lset,precision,[],mps,NsweepsMax);

[VZ,e,n]=expectationvalueTr(mps0,OsetZ);

[Vcurr,Tcur,]=expectationvalueTr(mps0,OsetCurrent);
[Vcurr2,Tcur2,n]=expectationvalueTr(mps0,OsetCurrent2);
current0=real((Tcur-Tcur2)/(N-1)/n);
Err0=imag((Tcur-Tcur2)/(N-1)/n);


fname = sprintf('DyVarJz%dW%dG%dDmax%dN%d.mat',Jz,Wd,Gamma,Dmax,N);
save(fname)
