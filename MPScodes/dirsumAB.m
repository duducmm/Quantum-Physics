function PSI=dirsumAB(PSI,PSI1)

    d1=size(PSI);

for j=1:d1(3);
PSI(:,:,j)=blkdiag(PSI(:,:,j),PSI1(:,:,j));
end



