function [n]=TrMPS(mps,N)
d=size(mps{1},3);
mpsid=eye(sqrt(d)); 
mpsid=reshape(mpsid,[1,1,d]); 
% trace
n=1;
X=eye(d); 
X=reshape(X,[1,1,d,d]);
for j=N:-1:1
n=updateCright(n,mpsid,X,mps{j}); 
end
