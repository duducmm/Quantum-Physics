function [E,mps]=minimizeE(hset,D,precision,mpsB)
[M,N]=size(hset); d=size(hset{1,1},1); 
mps=createrandommps(N,D,d); mps=prepare(mps);
% storage-initialization
Hstorage=initHstorage(mps,hset,d);
if ~isempty(mpsB), Cstorage=initCstorage(mps,[],mpsB,N); end 
P=[];
% optimization sweeps 
whileindex=0;
Evalues=[];

while 1
    clc
    whileindex=whileindex+1
std(Evalues)/abs(mean(Evalues))
Evalues=[];

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
[A,E]=minimizeE_onesite(hsetj,Hleft,Hright,P); 
[A,U]=prepare_onesite(A,'lr');
mps{j}=A;
Evalues=[Evalues,E];
% storage-update 
for m=1:M
h=reshape(hset{m,j},[1,1,d,d]);
Hstorage{m,j+1}=updateCleft(Hleft{m},A,h,A);
end
if ~isempty(mpsB) 
    Cstorage{j+1}=updateCleft(Cleft,A,[],B);
end
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
[A,E]=minimizeE_onesite(hsetj,Hleft,Hright,P); 
[A,U]=prepare_onesite(A,'rl');
mps{j}=A; Evalues=[Evalues,E];
% storage-update 
for m=1:M
h=reshape(hset{m,j},[1,1,d,d]);
Hstorage{m,j}=updateCright(Hright{m},A,h,A); 
end
if ~isempty(mpsB) 
    Cstorage{j}=updateCright(Cright,A,[],B);
end
end
if (std(Evalues)/abs(mean(Evalues))<precision || whileindex>10) 
    mps{1}=contracttensors(mps{1},3,2,U,2,1); 
    mps{1}=permute(mps{1},[1,3,2]);
break;

end
end
% ************************ one-site optimization **************************
function [A,E]=minimizeE_onesite(hsetj,Hleft,Hright,P)
DAl=size(Hleft{1},1); DAr=size(Hright{1},1); d=size(hsetj{1},1);
% calculation of Heff 
M=size(hsetj,1);
Heff=0; 
for m=1:M
Heffm=contracttensors(Hleft{m},3,2,Hright{m},3,2);
Heffm=contracttensors(Heffm,5,5,hsetj{m},3,3); 
Heffm=permute(Heffm,[1,3,5,2,4,6]); 
Heffm=reshape(Heffm,[DAl*DAr*d,DAl*DAr*d]); 
Heff=Heff+Heffm;
end
% projection on orthogonal subspace 
if ~isempty(P), 
    Heff=P'*Heff*P; end
% optimization
options.disp=0;  
[A,E]=eigs(Heff,1,'sr',options); 
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