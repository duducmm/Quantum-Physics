N=20;
D=20; precision=1e-10;
% XYZ Hamiltonian

Jz=1;
Wi=0;
hz=Wi*2*(rand(N,1)-0.5);
hz(1)=0.0001;
M=3*(N-1)+N;
hset=cell(M,N);
Oset=cell(N,N);
StagMag=cell(N,N);

sx=[0,1;1,0]; sy=[0,-1i;1i,0]; sz=[1,0;0,-1];
id=eye(2); 
for m=1:M, 
    for j=1:N, 
        hset{m,j}=id;               
    end
end
    for j=1:N, 
    for jj=1:N, 
        Oset{j,jj}=id;               
    end
    end
    
for j=1:N, 
for jj=1:N, 
        StagMag{j,jj}=id;               
end
end
for j=1:N
StagMag{j,j}=(-1)^j*sz/N;
end
   Stagmzvalues=[];
for jJz=1:50
    
    Jz=(jJz-1)/(50-1)*5;

for j=1:(N-1)
hset{3*(j-1)+1,j}=sx;
hset{3*(j-1)+1,j+1}=sx;
hset{3*(j-1)+2,j}=sy;
hset{3*(j-1)+2,j+1}=sy;
hset{3*(j-1)+3,j}=sqrt(Jz)*sz;
hset{3*(j-1)+3,j+1}=sqrt(Jz)*sz;
hset{3*(N-1)+j,j}=hz(j)*sz;

Oset{j,j}=sz;
end
hset{3*(N-1)+j+1,j+1}=hz(j+1)*sz;
Oset{j+1,j+1}=sz;

% ground state energy
randn('state',0) 
[E0,mps0]=minimizeE(hset,D,precision,[]); 
fprintf('E0 = %g\n',E0);
% % first excited state 
% [E1,mps1]=minimizeE(hset,D,precision,mps0); 
% fprintf('E1 = %g\n',E1);

[VEe,e,n]=expectationvalue2(mps0,Oset);
[VEe,Stagmz,n]=expectationvalue2(mps0,StagMag);
Stagmzvalues=[Stagmzvalues,Stagmz/n];
end


plot(Stagmzvalues)
