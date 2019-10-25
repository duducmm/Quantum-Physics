
clear
clc

N=20;
Dwarm=20;
Dmax=50;
precision=10^(-12);
Epsilon=10^(-8);
EpsilonMPO=10^(-8);

Jx=1;
Jy=1;
Jz=0.5;
h=0;
Wd=0;
Gamma=0.1;
P=1;
g=0.1;%JC couplin
k=0.5*g;
alpha=1;

dtwarm=0.01;
dt=0.001;

Ndtwarm=round(5/k/dtwarm);
Ndt=round(2/k/dt);
Fid=ones(1,Ndtwarm+Ndt);
% magnetization in z-direction
OsetZ=cell(N,N);
OsetS=cell(N,N);
VsetZ=cell(1,N);
OsetCurrent=cell(N-1,N);
d=2;
dl=d^2;
sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1];
s=[0 0 ;1 0];
Id=eye(dl);
id=eye(d);
for j=1:N,
    VsetZ{j}=reshape(sz,1,1,dl);
    for jj=1:N,
        OsetZ{j,jj}=Id;
        OsetS{j,jj}=Id;
    end
end
for j=1:N-1
    OsetZ{j,j}=kron(id,sz);
    OsetS{j,j}=kron(id,s);
end
OsetZ{N,N}=kron(id,sz);
OsetS{N,N}=kron(id,s);

% starting state
mps=cell(1,N);
for j=1:N
    state=[1; 0];
    state=state*state';
    %state=eye(d);
    state=state/(trace(state));
    state=reshape(state,dl,1);
    mps{j}=reshape(state,[1,1,dl]);
end
% %mps=createrandommps(N,Dmin,dl);
mps=prepare(mps);

hz= Wd.*(rand(N,1)-0.5)*2+h;

% time evolution
TT=[linspace(1,Ndtwarm,Ndtwarm)*dtwarm Ndtwarm*dtwarm+linspace(1,Ndt,Ndt)*dt];

for step=1:Ndtwarm
    clc
    step
    mpsi=mps;

    [~,S,n]=expectationvalueTr(mps,OsetS);
    S=S/n;
    alpha=alpha+(-1i*g*S-k/2*alpha+10^(-4))*dtwarm;

    [MPO]=buildMPO(N,dl,hz,Jx,Jy,Jz,Gamma,P,g,alpha,sx,sy,sz,s,Id,id,dtwarm,Dwarm,EpsilonMPO,precision);
    mps=MpoMps(MPO,mps,N);
    [mpsSVD]=mpsTruncSvd2(mps,Dwarm,Epsilon);
    mpsSVD=prepare(mpsSVD);
    [mps]=VarCompressMPS(mps,mpsSVD,precision);
        
    Fid(step)=abs(braket(mpsi,mps))^2;
    
    [VZ,e,n]=expectationvalueTr(mps,OsetZ);
    
    Zsurf(:,step)=VZ.'/n;
    Valpha(step)=alpha;

end

for step=1:Ndt
    
    clc
    step
    mpsi=mps;
    
    [~,S,n]=expectationvalueTr(mps,OsetS);
    S=S/n;
    alpha=alpha+(-1i*g*S-k/2*alpha+10^(-4))*dt;

    [MPO]=buildMPO(N,dl,hz,Jx,Jy,Jz,Gamma,P,g,alpha,sx,sy,sz,s,Id,id,dt,Dmax,EpsilonMPO,precision);
    
    mps=MpoMps(MPO,mps,N);
    [mpsSVD]=mpsTruncSvd2(mps,Dmax,Epsilon);
    mpsSVD=prepare(mpsSVD);
    [mps]=VarCompressMPS(mps,mpsSVD,precision);
       
    Fid(Ndtwarm+step)=abs(braket(mpsi,mps))^2;    
    
    [VZ,e,n]=expectationvalueTr(mps,OsetZ);
    
    Zsurf(:,Ndtwarm+step)=VZ.'/n;
    
end

for nn=1:N
    for j=1:N,
        osetg2{1,j}=Id;
    end
    osetg2{1,ceil(N/2)}=kron(id,sz);
    osetg2{1,nn}=kron(id,sz);
    if nn==ceil(N/2)
        osetg2{1,nn}=kron(id,sz^2);
    end
    [~,e,Tr]=expectationvalueTr(mps,osetg2);
    Spacialg2(nn)=e/Tr;
end
Spacialg2=Spacialg2'./(VZ.'/n*VZ(ceil(N/2)).'/n);

errorbar([0:1/(N-1):1],real(VZ.'/n),imag(VZ.'/n),'o-')
plot(Valpha.*conj(Valpha))
plot([1:1:N],real(Spacialg2),'o-')
plot([1:1:floor(N/2)],real(Spacialg2(ceil(N/2)+1:end)),'o-')

