function [n]=normMPS(mps,N)
d=size(mps{1},3);
% norm
n=1;
X=eye(d); 
X=reshape(X,[1,1,d,d]); 
for j=N:-1:1
n=updateCright(n,mps{j},X,mps{j}); 
end