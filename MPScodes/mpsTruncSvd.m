function [mps,K]=mpsTruncSvd(mps,Dmax,Epsilon)
N=length(mps);
for i=N:-1:2
    [mps{i},U]=DecimPrep(mps{i},'rl',Dmax,Epsilon);
    mps{i-1}=contracttensors(mps{i-1},3,2,U,2,1);
    mps{i-1}=permute(mps{i-1},[1,3,2]);
end
b=reshape(mps{1},[numel(mps{1}),1]);
K=b'*b;
mps{1}=mps{1}/sqrt(K);
