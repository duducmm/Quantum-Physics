function [E,mps]=minimizeHLL(hset,D,precision,mpsB,mps)
[M,N]=size(hset); d=size(hset{1,1},1); 
mpsR=createrandommps(N,D,d); 

for j=1:N
    d1=size(mpsR{j});
A=mps{j};
d2=size(mps{j});
A=padarray(A,d1-d2,0,'post');
mps{j}=A;
end
[n]=normMPS(mps,N);
for j=1:N
mps{j}=mps{j}/sqrt(abs(n))^(1/N);
end
mps=prepare(mps);
PSI=mps{1};
PSI=permute(PSI,[1,3,2]);

% storage-initialization
Hstorage=initLstorage(mps,hset,d);
if ~isempty(mpsB), Cstorage=initCstorage(mps,[],mpsB,N); end 
P=[];
% optimization sweeps 
whileindex=0;
Evalues=[];

while 1
    
    whileindex=whileindex+1
% std(Evalues)
% abs(mean(Evalues))
Evalues=[];
PSI=permute(PSI,[1,3,2]);

% ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
for j=1:(N-1)
% projector-calculation 
if ~isempty(mpsB)
B=mpsB{j};
Cleft=Cstorage{j};
Cright=Cstorage{j+1};
P=calcprojector_onesite(B,Cleft,Cright);
end
% optimization
Hleft=Hstorage(:,j);
Hright=Hstorage(:,j+1);
hsetj=hset(:,j);
PSI1=PSI;
[PSI,E]=minimizeL_onesite(hsetj,Hleft,Hright,P,PSI);
FID=abs(contracttensors(conj(PSI1),3,[1,2,3],PSI,3,[1,2,3]))

[A,U]=prepare_onesite(PSI,'lr');

mps{j}=A;
Evalues=[Evalues,abs(E)];
% storage-update 

for m=1:M
h=reshape(hset{m,j},[1,1,d,d]);
Hstorage{m,j+1}=updateCleft(Hleft{m},A,h,A);
end
if ~isempty(mpsB) 
    Cstorage{j+1}=updateCleft(Cleft,A,[],B);
end


PSI=contracttensors(conj(A),3,[1,3],PSI,3,[1,3]);
PSI=contracttensors(PSI,2,2,mps{j+1},3,1);


end
% ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
for j=N:(-1):2
% projector-calculation 
if ~isempty(mpsB)
B=mpsB{j};
Cleft=Cstorage{j};
Cright=Cstorage{j+1}; P=calcprojector_onesite(B,Cleft,Cright);
end
% minimization
Hleft=Hstorage(:,j);
Hright=Hstorage(:,j+1);
hsetj=hset(:,j); 
PSI=permute(PSI,[1,3,2]);
PSI1=PSI;
[PSI,E]=minimizeL_onesite(hsetj,Hleft,Hright,P,PSI);

FID=abs(contracttensors(conj(PSI1),3,[1,2,3],PSI,3,[1,2,3]))
[A,U]=prepare_onesite(PSI,'rl');

mps{j}=A; Evalues=[Evalues,abs(E)];
% storage-update 
for m=1:M
h=reshape(hset{m,j},[1,1,d,d]);
Hstorage{m,j}=updateCright(Hright{m},A,h,A); 
end
if ~isempty(mpsB) 
    Cstorage{j}=updateCright(Cright,mpsid,[],B);
end

 PSI=contracttensors(PSI,3,[2,3],conj(A),3,[2,3]);
 PSI=contracttensors(mps{j-1},3,2,PSI,2,1);

end
if ((std(Evalues)/sqrt(2*N))<precision || whileindex>4) 
    mps{1}=contracttensors(mps{1},3,2,U,2,1); 
    mps{1}=permute(mps{1},[1,3,2]);
break;

end
end
% ************************ one-site optimization **************************
function [A,E]=minimizeL_onesite(hsetj,Hleft,Hright,P,PSI)
DAl=size(Hleft{1},1); DAr=size(Hright{1},1); d=size(hsetj{1},1);
% calculation of Heff 
M=size(hsetj,1);
Heff=0; 
%tic
for m=1:M
Heffm=contracttensors(Hleft{m},3,2,Hright{m},3,2);
Heffm=contracttensors(Heffm,5,5,hsetj{m},3,3); 
Heffm=permute(Heffm,[1,3,5,2,4,6]); 
Heffm=reshape(Heffm,[DAl*DAr*d,DAl*DAr*d]); 
Heff=Heff+Heffm;
end
%toc
% projection on orthogonal subspace
if ~isempty(P), 
    Heff=P'*Heff*P; end
% optimization
options.disp=0; 
options.maxit=100;
options.tol=10^(3);
options.v0=reshape(PSI,[DAl*DAr*d,1]);
%tic
[A,E]=eigs(Heff,1,'sr',options); 
%toc
% V0q=zeros(size(Heff,1),1);
% V0q(DAl*DAr*d)=1;
% [A,E]= lsqr(Heff,V0q,10^(-14),10000);

% A=reshape(A,[sqrt(DAl*DAr*d),sqrt(DAl*DAr*d)]);
% A=(A'*A);
% A=reshape(A,[DAl*DAr*d,1]);
% A=A/norm(A);

%[A,E]=eigs(Heff,1,0,options); 
%A=exp(-1i*angle(A(1)))*A;
if ~isempty(P), 
    A=P*A; 
end
A=reshape(A,[DAl,DAr,d]);

function [P]=calcprojector_onesite(B,Cleft,Cright)
y=contracttensors(Cleft,3,3,B,3,1); 
y=contracttensors(y,4,[2,3],Cright,3,[2,3]);
y=permute(y,[1,3,2]);
y=reshape(y,[prod(size(y)),1]);
Q=orth([y,eye(size(y,1))]);
P=Q(:,2:end);