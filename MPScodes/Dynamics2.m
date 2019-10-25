clear 
clc
delete(gcp('nocreate'))

parpool('local',16)
 rng('shuffle')

N=20
D=2*N
precision=10^(-12); 
dt=0.001; 
Ndt=3*N/dt;
Ndt2=100/dt;
%Ndt*dt
%Ndt2*dt
Nrea=16;
Jx=1;
Jy=1;
Jz=1;
Wd=0;
hz1=0;
Gamma=1
% magnetization in z-direction
OsetCurrent=cell(N-1,N);
Oset=cell(N,N);
Jump1=cell(1,N);
Jump2=cell(1,N);
JumpP=cell(N,N);

sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1]; id=eye(2); 
s=[0 0 ; 1 0];
d=2;
% for j=1:N, 
%     oset{1,j}=id; 
% end
% oset{1,jflipped}=sz;

 for j=1:N,
     Jump1{1,j}=reshape(id,[1,1,2,2]);     
     Jump2{1,j}=reshape(id,[1,1,2,2]);
    for jj=1:N, 
        Oset{j,jj}=id; 
        OsetCurrent{j,jj}=id;
        JumpP{j,jj}=id;
    end   
 end
 OsetCurrent2=OsetCurrent;
 jreno=10^0;
 Jump1{1,1}=jreno*reshape(s',[1,1,2,2]);
 Jump2{1,N}=jreno*reshape(s,[1,1,2,2]);
 JumpP{1,1}=Gamma*dt*(s*s');
 JumpP{N,N}=Gamma*dt*(s'*s);
 
for j=1:N-1
Oset{j,j}=sz;
OsetCurrent{j,j}=sx;
OsetCurrent{j,j+1}=sy;
OsetCurrent2{j,j}=sy;
OsetCurrent2{j,j+1}=sx;

end
Oset{N,N}=sz;

% time evolution operator
I=reshape(id,[1,1,2,2]);

Mzmed=zeros(Nrea,N);
SCurrent=zeros(Nrea,1);
parfor rea=1:Nrea
mpo_even=cell(1,N);
mpo_odd=cell(1,N);
mpo_local=cell(1,N);

hz= Wd.*(rand(N,1)-0.5)*2+hz1;

    for j=1:N
    mpo_even{j}=I; 
    mpo_odd{j}=I; 
    h=hz(j)*sz;
    mpo_local{j}=reshape(expm(-1i*dt/2*h),[1,1,d,d]);
    end
       h=hz(1)*sz-1i*Gamma*(s*s')/2;
    mpo_local{1}=reshape(expm(-1i*dt/2*h),[1,1,d,d]);
       h=hz(N)*sz-1i*Gamma*(s'*s)/2;
    mpo_local{N}=reshape(expm(-1i*dt/2*h),[1,1,d,d]);
    
    
    
for j=1:2:(N-1)
  
h=Jx*kron(sx,sx)+Jy*kron(sy,sy)+Jz*kron(sz,sz);
w=expm(-1i*dt*h);
w=reshape(w,[2,2,2,2]); 
w=permute(w,[1,3,2,4]);
w=reshape(w,[4,4]); [U,S,V]=svd2(w); eta=size(S,1);
U=U*sqrt(S); V=sqrt(S)*V;
U=reshape(U,[2,2,eta]);
U=permute(U,[4,3,2,1]); 
V=reshape(V,[eta,2,2]);
V=permute(V,[1,4,3,2]);

    mpo_odd{j}=U;
    mpo_odd{j+1}=V;
end 
for j=2:2:(N-1)
    
h=Jx*kron(sx,sx)+Jy*kron(sy,sy)+Jz*kron(sz,sz);
w=expm(-1i*dt*h/2);
w=reshape(w,[2,2,2,2]); 
w=permute(w,[1,3,2,4]);
w=reshape(w,[4,4]); [U,S,V]=svd2(w); eta=size(S,1);
U=U*sqrt(S); V=sqrt(S)*V;
U=reshape(U,[2,2,eta]);
U=permute(U,[4,3,2,1]); 
V=reshape(V,[eta,2,2]);
V=permute(V,[1,4,3,2]);

    mpo_even{j}=U;
    mpo_even{j+1}=V;
end
% starting state 
[mps]=createrandommps(N,D,d);
[VEe,mz,n]=expectationvalue2(mps,Oset);
mps{1}=mps{1}/sqrt(abs(n));

% mps0=cell(1,N);
% for j=1:N
%     state=[0; 1];
% mps0{j}=reshape(state,[1,1,2]);
% 
% end
% for j=2:2:N
%     state=[1; 0];
%     mps0{j}=reshape(state,[1,1,2]);
% 
% end

% time evolution
mzvalues=zeros(Ndt2,N);
cur=zeros(Ndt2,1);
mps2=mps;

TT=linspace(1,Ndt,Ndt)*dt;
for step1=1:Ndt
  
[mps,]=reduceD(mps2,mpo_local,D,precision); 
[mps,]=reduceD(mps,mpo_even,D,precision); 
[mps,]=reduceD(mps,mpo_odd,D,precision);
[mps,]=reduceD(mps,mpo_even,D,precision); 
[mps,]=reduceD(mps,mpo_local,D,precision); 
[VEe,~,n]=expectationvalue2(mps,JumpP);
p1=max(real(VEe(1)/n),0);
p2=max(real(VEe(N)/n),0);

Oracle=rand;
if Oracle<p1
[mps,]=reduceD(mps,Jump1,D,precision); 
elseif (p1<Oracle)&&(Oracle<(p1+p2))
[mps,]=reduceD(mps,Jump2,D,precision); 
end

[VEe,~,n]=expectationvalue2(mps,JumpP);
p1=max(real(VEe(1)/n),0);
p2=max(real(VEe(N)/n),0);

Oracle=rand;
if Oracle<p1
[mps]=JumpsMps(mps,Jump1,1);
elseif (p1<Oracle)&&(Oracle<(p1+p2))
[mps]=JumpsMps(mps,Jump2,N);
end

mps2=mps;

end
for step=1:Ndt2
  
[mps,]=reduceD(mps2,mpo_local,D,precision); 
[mps,]=reduceD(mps,mpo_even,D,precision); 
[mps,]=reduceD(mps,mpo_odd,D,precision);
[mps,]=reduceD(mps,mpo_even,D,precision); 
[mps,]=reduceD(mps,mpo_local,D,precision); 
[VEe,~,n]=expectationvalue2(mps,JumpP);
p1=max(real(VEe(1)/n),0);
p2=max(real(VEe(N)/n),0);

Oracle=rand;
if Oracle<p1
[mps,]=reduceD(mps,Jump1,D,precision); 
elseif (p1<Oracle)&&(Oracle<(p1+p2))
[mps,]=reduceD(mps,Jump2,D,precision); 
end

mps2=mps;

[VEe,mz,]=expectationvalue2(mps2,Oset);
[Vcurr,Tcur,]=expectationvalue2(mps2,OsetCurrent);
[Vcurr2,Tcur2,n]=expectationvalue2(mps2,OsetCurrent2);
Vcurr=Vcurr-Vcurr2;
cur(step)=(Tcur-Tcur2)/(N-1);
mzvalues(step,:)=real(VEe);

end


Mzmed(rea,:)=mean(mzvalues,1);
SCurrent(rea)=mean(cur);
end

Prof=mean(Mzmed,1);
currentR=Gamma*(Prof(N)+1)/2;
current=mean(real(SCurrent),1);
Err=std(real(SCurrent),0,1)/sqrt(Nrea);
ErrR=std(Mzmed,0,1)/sqrt(Nrea);


fname = sprintf('MPSJz%dW%dG%dD%dN%d.mat',Jz,Wd,Gamma,D,N);
save(fname)
delete(local)

