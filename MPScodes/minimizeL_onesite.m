function [A,E]=minimizeL_onesite(Lsetj,Hleft,Hright,P,PSI)
DAl=size(Hleft{1},1); 
DAr=size(Hright{1},1);
d=size(Lsetj{1},1);
% calculation of Heff
M=size(Lsetj,1);
Heff=0; 
%tic
for m=1:M
Heffm=contracttensors(Hleft{m},3,2,Hright{m},3,2);
Heffm=contracttensors(Heffm,5,5,Lsetj{m},3,3); 
Heffm=permute(Heffm,[1,3,5,2,4,6]); 
Heffm=reshape(Heffm,[DAl*DAr*d,DAl*DAr*d]);
Heff=Heff+Heffm;
end
%toc
% projection on orthogonal subspace
if ~isempty(P), 
    Heff=P'*Heff*P; end
% optimization
% options.disp=0;
% options.isreal=0;
% options.maxit=300;
 options.tol=10^(-17);
%options.p=40;
options.v0=reshape(PSI,[DAl*DAr*d,1]);
%tic
[A,E]=eigs(Heff,1,0,options); 
%toc
if ~isempty(P), 
    A=P*A; 
end
A=reshape(A,[DAl,DAr,d]);