function [Hstorage]=initHstorage(mps,hset,d)
[M,N]=size(hset);
Hstorage=cell(M,N+1);
for m=1:M,
    Hstorage{m,1}=1; 
    Hstorage{m,N+1}=1;
end
for j=N:-1:2
for m=1:M
h=reshape(hset{m,j},[1,1,d,d]); 
Hstorage{m,j}=updateCright(Hstorage{m,j+1},mps{j},h,mps{j});
end
end


function [Cright]=updateCright(Cright,B,X,A)
if isempty(X), 
    X=reshape(eye(size(B,3)),[1,1,2,2]); 
end
Cright=contracttensors(A,3,2,Cright,3,3); 
Cright=contracttensors(X,4,[2,4],Cright,4,[4,2]);
Cright=contracttensors(conj(B),3,[2,3],Cright,4,[4,2]);
