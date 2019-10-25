function [Cstorage]=initCstorage2(mpsB,X,mpsA,N)
Cstorage=cell(1,N+1); Cstorage{1}=1; Cstorage{N+1}=1;
for i=N:-1:2
Cstorage{i}=updateCright(Cstorage{i+1},mpsB{i},X,mpsA{i}); 
end