function [F]=Fidelity(mps,mps0)
d=size(mps{1},3);
[~,N]=size(mps); 

% norm
F=1;
X=eye(d); 
X=reshape(X,[1,1,d,d]); 
for j=N:-1:1
F=updateCright(F,mps0{j},X,mps{j}); 
end
