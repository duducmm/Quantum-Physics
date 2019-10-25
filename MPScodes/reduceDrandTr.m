function [mpsB,n]=reduceDrandTr(mpsA,mpoX,DB,precision)
N=length(mpsA);
d=size(mpsA{1},3); 
mpsB=createrandommps(N,DB,d);
mpsB=prepare(mpsB);
n=TrMPS(mpsB,N);
mpsB{1}=mpsB{1}/n;
% initialization of the storage 
Cstorage=initCstorage(mpsB,mpoX,mpsA,N);
% optimization sweeps 
while 1
Kvalues=[];
% ****************** cycle 1: j -> j+1 (from 1 to N-1) **************** 
for j=1:(N-1)
% optimization
Cleft=Cstorage{j};
Cright=Cstorage{j+1};
A=mpsA{j}; X=mpoX{j};
[B,K]=reduceD2_onesite(A,X,Cleft,Cright);
[B,U]=prepare_onesite(B,'lr');
mpsB{j}=B;
Kvalues=[Kvalues,K];
% storage-update
Cstorage{j+1}=updateCleft(Cleft,B,X,A);
end
% ****************** cycle 2: j -> j-1 (from N to 2) ****************** 
for j=N:(-1):2
% optimization
Cleft=Cstorage{j};
Cright=Cstorage{j+1};
A=mpsA{j}; X=mpoX{j}; 
[B,K]=reduceD2_onesite(A,X,Cleft,Cright);
[B,U]=prepare_onesite(B,'rl');
mpsB{j}=B;
Kvalues=[Kvalues,K];
% storage-update
Cstorage{j}=updateCright(Cright,B,X,A);
end
if std(Kvalues)/mean(Kvalues)<precision
    mpsB{1}=contracttensors(mpsB{1},3,2,U,2,1);
    n=TrMPS(mpsB,N);
break
end
end


% ************************ one-site optimization ************************** 
function [B,K]=reduceD2_onesite(A,X,Cleft,Cright)
Cleft=contracttensors(Cleft,3,3,A,3,1);
Cleft=contracttensors(Cleft,4,[2,4],X,4,[1,4]);
B=contracttensors(Cleft,4,[3,2],Cright,3,[2,3]);
B=permute(B,[1,3,2]);
b=reshape(B,[numel(B),1]);
K=b'*b;



