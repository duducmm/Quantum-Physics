parpool('local',16)

rng('shuffle')

N=20;

Jx=1;
Jy=1;
Jz=1;
h=0;
Wd=0;
Gamma=1;
mu=0.1;

Dwarm=20;
Dmax=100;
precision=10^(-8);
Epsilon=10^(-8);
precisionMPO=10^(-12);
EpsilonMPO=10^(-8);
Sweepmax=10;
SweepmaxMPO=100;

dtwarm=0.1;
dt=0.05;
Ndtwarm=round(2*N/dtwarm);
Ndt=round(4*N/dt);

Nrea=16;

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

Zvec=zeros(Nrea,N);
CellCurrent=zeros(1,Nrea);
CellErr=zeros(1,Nrea);

parfor rea=1:Nrea
    realizacao=rea
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
    mps=CanonizeA(mps);
    
    hz= Wd.*(rand(N,1)-0.5)*2+h;
    
    
    % time evolution operator
    dt1=dtwarm/(4-4^(1/3));
    dt2=dt1;
    dt3=dtwarm-2*dt1-2*dt2;
    [MPO1]=UMPOXXZbd(N,dl,hz,Jx,Jy,Jz,Gamma,mu,sx,sy,sz,s,Id,id,dt1,Dwarm,EpsilonMPO,precisionMPO,SweepmaxMPO);
    [MPO2]=UMPOXXZbd(N,dl,hz,Jx,Jy,Jz,Gamma,mu,sx,sy,sz,s,Id,id,dt2,Dwarm,EpsilonMPO,precisionMPO,SweepmaxMPO);
    [MPO3]=UMPOXXZbd(N,dl,hz,Jx,Jy,Jz,Gamma,mu,sx,sy,sz,s,Id,id,dt3,Dwarm,EpsilonMPO,precisionMPO,SweepmaxMPO);
    
    % time evolution
    TT=[linspace(1,Ndtwarm,Ndtwarm)*dtwarm Ndtwarm*dtwarm+linspace(1,Ndt,Ndt)*dt];
    
    for step=1:Ndtwarm
        mps=TrotterSuzuki4th(MPO1,MPO2,MPO3,mps,N,Dwarm,Epsilon,precision,Sweepmax);
    end
    
    dt1=dt/(4-4^(1/3));
    dt2=dt1;
    dt3=dt-2*dt1-2*dt2;
    [MPO1]=UMPOXXZbd(N,dl,hz,Jx,Jy,Jz,Gamma,mu,sx,sy,sz,s,Id,id,dt1,Dmax,EpsilonMPO,precisionMPO,SweepmaxMPO);
    [MPO2]=UMPOXXZbd(N,dl,hz,Jx,Jy,Jz,Gamma,mu,sx,sy,sz,s,Id,id,dt2,Dmax,EpsilonMPO,precisionMPO,SweepmaxMPO);
    [MPO3]=UMPOXXZbd(N,dl,hz,Jx,Jy,Jz,Gamma,mu,sx,sy,sz,s,Id,id,dt3,Dmax,EpsilonMPO,precisionMPO,SweepmaxMPO);
    
    for step=1:Ndt
        mps=TrotterSuzuki4th(MPO1,MPO2,MPO3,mps,N,Dmax,Epsilon,precision,Sweepmax);
    end
    
    [VZ,e,n]=expectationvalueTr(mps,OsetZ);
    
    Zvec(rea,:)=VZ.'/n;
    
    [Vcurr,Tcur,]=expectationvalueTr(mps,OsetCurrent);
    [Vcurr2,Tcur2,n]=expectationvalueTr(mps,OsetCurrent2);
    current=real((Tcur-Tcur2)/(N-1)/n);
    Err=imag((Tcur-Tcur2)/(N-1)/n);
    CellCurrent(rea)=current;
    CellErr(rea)=Err;
    
end

current=mean(CellCurrent);
Err=mean(CellErr);
Perfil=mean(Zvec,1);
ErrPerfil=std(real(Zvec),0,1)/sqrt(Nrea);
ErrCurrent=std(real(CellCurrent),0,2)/sqrt(Nrea);

fname = sprintf('DMdyJz%dW%dG%dmu%dDmax%dN%d.mat',Jz,Wd,Gamma,mu,Dmax,N);
save(fname)
