function [mpsB]=mpsCompressSvd(mps,Epsilon,Dmax,N)

[mpsB{N},US] = canonizeA(mps{N},Dmax,Epsilon);
for site = (N-1):(-1):1
    mpsB{site} = permute(contract(mps{site},2,US,1),[1 3 2]);
    [mpsB{site},US] = canonizeA(mps{site},Dmax,Epsilon);
end
mpsB{1} = mpsB{1}*sign(US);