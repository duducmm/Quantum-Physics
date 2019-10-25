clear
clc

N=20;
D=20;
precision=10^(-8);
NsweepsMax=8;
% XYZ Hamiltonian
Gamma=1;
GAMMA=zeros(N,1);
GAMMA(1)=Gamma;
GAMMA(N)=Gamma;

Jx=1;
Jy=Jx;
Jz=10;
Jx2=0;

h=0;
Wi=0;
hz=Wi*2*(rand(N,1)-0.5)+h;
M=N+3*(N-2)+3+3*(N-2)+3;%+2*(N-2);
%M=N+1*(N-2)+1+1*(N-2)+1;
%M=N;
Lset=cell(M,N);
oset=cell(1,N);
J=[0 0; 1 0]';

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

for m=1:M, 
    for j=1:N, 
        Lset{m,j}=Id;               
    end
end

    hz(:)=h;
    hz(1)=hz(1)+0.0001;

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
mps=createrandommps(N,D,dl); 
mps=prepare(mps);

[E0,mps]=minimizeL(Lset,precision,[],mps,NsweepsMax);


[VZ,e,~]=expectationvalueTr(mps,OsetZ);

[Vcurr,Tcur,]=expectationvalueTr(mps,OsetCurrent);
[Vcurr2,Tcur2,n]=expectationvalueTr(mps,OsetCurrent2);
current=real((Tcur-Tcur2)/(N-1)/n);
Err=imag((Tcur-Tcur2)/(N-1)/n);
