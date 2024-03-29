function [VEe,e,n]=expectationvalueTr(mps,hset)
[M,N]=size(hset); 
d=size(mps{1},3);
% expectation value 
mpsid=eye(sqrt(d)); 
mpsid=reshape(mpsid,[1,1,d]); 
e=0;
for m=1:M
em=1;
for j=N:-1:1
h=hset{m,j};
h=reshape(h,[1,1,d,d]); em=updateCright(em,mpsid,h,mps{j});
end
e=e+em;
VEe(m)=em;
end
% trace
n=1;
X=eye(d); 
X=reshape(X,[1,1,d,d]);

for j=N:-1:1
n=updateCright(n,mpsid,X,mps{j}); 
end
%e=e/n;