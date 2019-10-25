function [mps,K]=prepare(mps)
N=length(mps);
for i=N:-1:2
    [mps{i},U]=prepare_onesite(mps{i},'rl');
    mps{i-1}=contracttensors(mps{i-1},3,2,U,2,1);
    mps{i-1}=permute(mps{i-1},[1,3,2]);
end
b=reshape(mps{1},[numel(mps{1}),1]);
K=b'*b;
mps{1}=mps{1}/sqrt(K);
function [B,U,DB]=prepare_onesite(A,direction)
[D1,D2,d]=size(A);
switch direction
    case 'lr'
        A=permute(A,[3,1,2]); A=reshape(A,[d*D1,D2]); [B,S,U]=svd2(A); DB=size(S,1);
        B=reshape(B,[d,D1,DB]); B=permute(B,[2,3,1]);
        U=S*U;
    case 'rl'
        A=permute(A,[1,3,2]); A=reshape(A,[D1,d*D2]); [U,S,B]=svd2(A); DB=size(S,1);
        B=reshape(B,[DB,d,D2]); B=permute(B,[1,3,2]); 
        U=U*S;
end
