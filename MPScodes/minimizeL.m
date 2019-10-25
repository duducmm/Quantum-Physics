function [E,mps]=minimizeL(Lset,precision,mpsB,mps,NsweepsMax)
[M,N]=size(Lset); d=size(Lset{1,1},1); 

PSI=mps{1};
PSI=permute(PSI,[1,3,2]);

% storage-initialization
Hstorage=initLstorage(mps,Lset,d);
if ~isempty(mpsB), Cstorage=initCstorage(mps,[],mpsB,N); end 
P=[];
% optimization sweeps 
Sweep=0;
 
while 1
    Evalues=[];

        Sweep=Sweep+1
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
Lsetj=Lset(:,j);
PSI1=PSI;
[PSI,E]=minimizeL_onesite(Lsetj,Hleft,Hright,P,PSI);
FID=1-abs(contracttensors(conj(PSI1),3,[1,2,3],PSI,3,[1,2,3]))

[A,U]=prepare_onesite(PSI,'lr');
mps{j}=A;
Evalues=[Evalues,abs(E)];
% storage-update 

for m=1:M
h=reshape(Lset{m,j},[1,1,d,d]);
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
Lsetj=Lset(:,j); 
PSI=permute(PSI,[1,3,2]);
PSI1=PSI;

[PSI,E]=minimizeL_onesite(Lsetj,Hleft,Hright,P,PSI);
FID=1-abs(contracttensors(conj(PSI1),3,[1,2,3],PSI,3,[1,2,3]))

[A,U]=prepare_onesite(PSI,'rl');

mps{j}=A; Evalues=[Evalues,abs(E)];
% storage-update 
for m=1:M
h=reshape(Lset{m,j},[1,1,d,d]);
Hstorage{m,j}=updateCright(Hright{m},A,h,A); 
end
if ~isempty(mpsB) 
    Cstorage{j}=updateCright(Cright,mpsid,[],B);
end

 PSI=contracttensors(PSI,3,[2,3],conj(A),3,[2,3]);
 PSI=contracttensors(mps{j-1},3,2,PSI,2,1);

end
if (std(Evalues)/mean(Evalues)<precision || Sweep>=NsweepsMax) 
mps{1}=contracttensors(mps{1},3,2,U,2,1); 
    mps{1}=permute(mps{1},[1,3,2]);
    b=reshape(mps{1},[numel(mps{1}),1]);
K=b'*b;
mps{1}=mps{1}/sqrt(K);
break;

end

end
% ************************ one-site optimization **************************
function [A,E]=minimizeL_onesite(Lsetj,Hleft,Hright,P,PSI)
DAl=size(Hleft{1},1); 
DAr=size(Hright{1},1);
d=size(Lsetj{1},1);
% calculation of Heff
M=size(Lsetj,1);
Heff=0; 
%tic
for m=1:M
Heffm=contracttensors(Hleft{m},3,2,Hright{m},3,2);
Heffm=contracttensors(Heffm,5,5,Lsetj{m},3,3); 
Heffm=permute(Heffm,[1,3,5,2,4,6]); 
Heffm=reshape(Heffm,[DAl*DAr*d,DAl*DAr*d]);
Heff=Heff+Heffm;
end
%toc
% projection on orthogonal subspace
if ~isempty(P), 
    Heff=P'*Heff*P; end
% optimization
% options.disp=0;
% options.isreal=0;
% options.maxit=300;
 options.tol=10^(-17);
%options.p=40;
options.v0=reshape(PSI,[DAl*DAr*d,1]);
%tic
[A,E]=eigs(Heff,1,0,options); 
%toc
if ~isempty(P), 
    A=P*A; 
end
A=reshape(A,[DAl,DAr,d]);

function [P]=calcprojector_onesite(B,Cleft,Cright)
y=contracttensors(Cleft,3,3,B,3,1); 
y=contracttensors(y,4,[2,3],Cright,3,[2,3]);
y=permute(y,[1,3,2]);
y=reshape(y,[numel(y),1]);
Q=orth([y,eye(size(y,1))]);
P=Q(:,2:end);