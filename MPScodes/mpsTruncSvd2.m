function [mps,K]=mpsTruncSvd2(mps,Dmax,Epsilon)
N=length(mps);
for i=1:1:N-1
    [mps{i},U]=DecimPrep(mps{i},'lr',Dmax,Epsilon);
    mps{i+1}=contracttensors(U,2,2,mps{i+1},3,1);
end
b=reshape(mps{N},[numel(mps{N}),1]);
K=b'*b;
mps{N}=mps{N}/sqrt(K);