clear 
clc


N=20
D=max(2*N,50)
precision=1e-12; 
Jx=1;
Jy=1;
Jz=1;
h=0;
Wd=0;
Gamma=1;
dt=0.01; 
Ndt=max(4*N/dt,100/dt);
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
% state=[1; 0];
 %state=state*state';
 state=eye(d); 
 state=state/(trace(state));
 state=reshape(state,dl,1);
mps{j}=reshape(state,[1,1,dl]);
end


%end
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

    w=expm(dt*LH);
    
    if j==1
    w=expm(dt*(LH+H1+L1));
    elseif j==(N-1)
        w=expm(dt*(LH+HN+LN));
    end
    
w=reshape(w,[dl,dl,dl,dl]);
w=permute(w,[4,2,3,1]); 
w=reshape(w,[dl^2,dl^2]); 
[U,S,V]=svd2(w); 
eta=size(S,1);
U=U*sqrt(S); 
V=sqrt(S)*V;
U=reshape(U,[dl,dl,eta]);
U=permute(U,[4,3,2,1]);
V=reshape(V,[eta,dl,dl]);
V=permute(V,[1,4,3,2]);

    

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
w=reshape(w,[dl,dl,dl,dl]);
w=permute(w,[4,2,3,1]); 
w=reshape(w,[dl^2,dl^2]); 
[U,S,V]=svd2(w); 
eta=size(S,1);
U=U*sqrt(S); 
V=sqrt(S)*V;
U=reshape(U,[dl,dl,eta]);
U=permute(U,[4,3,2,1]);
V=reshape(V,[eta,dl,dl]);
V=permute(V,[1,4,3,2]);
    mpo_even{j}=U;
    mpo_even{j+1}=V;
end


% time evolution
TT=linspace(1,Ndt,Ndt)*dt;
for step=1:Ndt
   [mps]=reduceD(mps,mpo_even,D,precision); 
    [mps]=reduceD(mps,mpo_odd,D,precision);
    [mps,K]=reduceD(mps,mpo_even,D,precision); 
    clc
    step


end
[VZ,e,n]=expectationvalueTr(mps,OsetZ);

Zsurf(:,step)=VZ.'/n;

[Vcurr,Tcur,]=expectationvalueTr(mps,OsetCurrent);
[Vcurr2,Tcur2,n]=expectationvalueTr(mps,OsetCurrent2);
current=real((Tcur-Tcur2)/(N-1)/n);
Err=imag((Tcur-Tcur2)/(N-1)/n);

fname = sprintf('MPOJz%dW%dG%dD%dN%d.mat',Jz,Wd,Gamma,D,N);
save(fname)

