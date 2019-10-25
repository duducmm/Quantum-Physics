function [B,U,DB]=Canonize_onesite(A,direction)
[D1,D2,d]=size(A);
switch direction
    case 'lr'
        A=permute(A,[3,1,2]); A=reshape(A,[d*D1,D2]); [B,U]=qr2(A); DB=size(B,2);
        B=reshape(B,[d,D1,DB]); B=permute(B,[2,3,1]);
    case 'rl'
        A=permute(A,[1,3,2]); A=reshape(A,[D1,d*D2]); [B,U]=qr2(A'); B=B';U=U'; DB=size(B,1);
        B=reshape(B,[DB,d,D2]); B=permute(B,[1,3,2]);
end