clear
clc
warning('off','all')

N=50;
Dmax=60;
Dmin=5;
D=Dmin; precision=1e-4;
NsweepsMax=5;
% XYZ Hamiltonian
gridE=61;
Gamma=0.1;
Jz=0;
Jx=1;
Jx2=0.5;
Jy=0;
h=-3;
Wi=0;
hz=Wi*2*(rand(N,1)-0.5)+h;
hz(1)=0.0001;
%M=N+3*(N-2)+3+3*(N-2)+3%+2*(N-2);
M=N+1*(N-2)+1+1*(N-2)+1;
%M=N;
Lset=cell(M,N);
oset=cell(1,N);

sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1];
J=sqrt(Gamma)*[0 0; 1 0];
id=eye(2^2); 

Vec0=zeros(2^2,1);
Vec0(1)=1;

Vec0=id(:,2^2);


for m=1:M, 
    for j=1:N, 
        Lset{m,j}=id;               
    end
end

%  HLL=cell(M^2,N);
% for m=1:M^2, 
%     for j=1:N, 
%         HLL{m,j}=id;               
%     end
% end

%     for j=1:N, 
%     for jj=1:N, 
%         Oset{j,jj}=id;               
%     end
%     end
     for j=1:N, 
         oset{1,j}=id;               
     end
       
id=eye(2); 





for jh=1:1
    
    h=(-3+6*(jh-1)/(gridE-1));
    hz(:)=1;
    hz(1)=hz(1)+0.0001;
    Vjh(jh)=h;

for j=1:(N-1)
% Lset{N+3*(j-1)+1,j}=-1i*sqrt(Jx)*kron(id,sx);
% Lset{N+3*(j-1)+1,j+1}=sqrt(Jx)*kron(id,sx);
% Lset{N+3*(j-1)+2,j}=-1i*sqrt(Jy)*kron(id,sy);
% Lset{N+3*(j-1)+2,j+1}=sqrt(Jy)*kron(id,sy);
% Lset{N+3*(j-1)+3,j}=-1i*sqrt(Jz)*kron(id,sz);
% Lset{N+3*(j-1)+3,j+1}=sqrt(Jz)*kron(id,sz);
% 
% Lset{N+3*(N-2)+3+3*(j-1)+1,j}=1i*sqrt(Jx)*kron(sx.',id);
% Lset{N+3*(N-2)+3+3*(j-1)+1,j+1}=sqrt(Jx)*kron(sx.',id);
% Lset{N+3*(N-2)+3+3*(j-1)+2,j}=1i*sqrt(Jy)*kron(sy.',id);
% Lset{N+3*(N-2)+3+3*(j-1)+2,j+1}=sqrt(Jy)*kron(sy.',id);
% Lset{N+3*(N-2)+3+3*(j-1)+3,j}=1i*sqrt(Jz)*kron(sz.',id);
% Lset{N+3*(N-2)+3+3*(j-1)+3,j+1}=sqrt(Jz)*kron(sz.',id);

Lset{N+1*(j-1)+1,j}=-1i*sqrt(Jx)*kron(id,sx);
Lset{N+1*(j-1)+1,j+1}=sqrt(Jx)*kron(id,sx);
Lset{N+1*(N-2)+1+1*(j-1)+1,j}=1i*sqrt(Jx)*kron(sx.',id);
Lset{N+1*(N-2)+1+1*(j-1)+1,j+1}=sqrt(Jx)*kron(sx.',id);

Lset{j,j}=kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id))...
    -1i*hz(j)*kron(id,sz)+1i*hz(j)*kron(sz.',id);

%Lset{N+1,j}=Vec0*reshape(id,1,2^2);

end
j=j+1;
 Lset{j,j}=kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id))...
     -1i*hz(j)*kron(id,sz)+1i*hz(j)*kron(sz.',id);

% for j=1:(N-2)
% Lset{N+3*(N-2)+3+3*(N-2)+3+(j-1)+1,j}=-1i*sqrt(Jx2)*kron(id,sx);
% Lset{N+3*(N-2)+3+3*(N-2)+3+(j-1)+1,j+2}=sqrt(Jx2)*kron(id,sx);
% Lset{N+3*(N-2)+3+3*(N-2)+3+N-2+(j-1)+1,j}=1i*sqrt(Jx2)*kron(sx.',id);
% Lset{N+3*(N-2)+3+3*(N-2)+3+N-2+(j-1)+1,j+2}=sqrt(Jx2)*kron(sx.',id);
% end
 
 %Lset{N+1,j}=Vec0*reshape(id,1,2^2);
% for m1=1:M
%     for m2=1:M
%     for j=1:N
%        HLL{(m1-1)*M+m2,j}=Lset{m1,j}'*Lset{m2,j};
%     end
%     end
% end

% Asymptotic State
D=Dmin;
d=size(oset{1,1},1); 
mps0=createrandommps(N,D,d); 

for j=1:N
    d1=size(mps0{j});
 state=[1; 1];
 state=state*state';
 state=eye(sqrt(d)); 
 state=state/(trace(state));
 state=reshape(state,d,1);
mpsid1=reshape(state,[1,1,d]);
mpsid=zeros(d1);
for jj=1:d1(3)
mpsid(:,:,jj)=mpsid1(1,1,jj)*1;%eye(d1(1),d1(2));
end
mps0{j}=mpsid;
end
[n]=normMPS(mps0,N);
for j=1:N
mps0{j}=mps0{j}/sqrt(abs(n))^(1/N);
end
mps0=prepare(mps0);



while 1
[E0,mps0]=minimizeL(Lset,D,Dmin,precision/100,[],mps0,...
    NsweepsMax,hz,Gamma,Jx,Jy,Jz);
%[E0,mps0]=minimizeHLL(HLL,D,precision/10,[],mps0); 

if (abs(E0)<=precision || D>=Dmax)
%    if abs(E0)<10^(-4)
%        D=D-5;
%    end
    break
else 
    D=D+max(round(D/2),10);   
if D>Dmax
    D=Dmax;
end

end
end

if D<Dmin
    D=Dmin;
end


for nn=1:N/2
for j=1:N, 
         oset{1,j}=eye(4);                   
 end
oset{1,round(N/2)}=kron(eye(2),sx);
oset{1,round(N/2)+nn}=kron(eye(2),sx);
[~,e,Tr]=expectationvalueTr(mps0,oset);
SpacialXX(nn)=e/Tr;
end

clc
jh
E0
e/Tr
D


save('VarMPO50SitesGamma0p1h1')
end

