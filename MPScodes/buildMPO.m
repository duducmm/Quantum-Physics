function [MPO]=buildMPO(N,dl,hz,Jx,Jy,Jz,Gamma,P,g,alpha,sx,sy,sz,s,Id,id,dt,Dmax,EpsilonMPO,precision)

Dij=-1i*g*kron(kron(id,alpha'*s+alpha*s'),Id)/2+1i*g*kron(kron(alpha'*s+alpha*s',id),Id)/2 ...
    -1i*g*kron(Id,kron(id,alpha'*s+alpha*s'))/2+1i*g*kron(Id,kron(alpha'*s+alpha*s',id))/2;

H1=-1i*hz(1)*kron(kron(id,sz),Id)/2+1i*hz(1)*kron(kron(sz,id),Id)/2 ...
    -1i*g*kron(kron(id,alpha'*s+alpha*s'),Id)/2+1i*g*kron(kron(alpha'*s+alpha*s',id),Id)/2;
%g*(alpha'*s+alpha*s')
HN=-1i*hz(N)*kron(Id,kron(id,sz))/2+1i*hz(N)*kron(Id,kron(sz,id))/2 ...
    -1i*g*kron(Id,kron(id,alpha'*s+alpha*s'))/2+1i*g*kron(Id,kron(alpha'*s+alpha*s',id))/2;

J=sqrt(Gamma)*s;
L1=kron(kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)),Id);
J=sqrt(P)*s';
L1=L1+kron(kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)),Id);
L1=L1/2;

J=sqrt(Gamma)*s;
LN=kron(Id,kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)));
J=sqrt(P)*s';
LN=LN+kron(Id,kron(conj(J),J)-1/2*(kron(id,J'*J)+kron((J'*J).',id)));
LN=LN/2;

Lij=L1+LN+Dij;



I=reshape(Id,[1,1,dl,dl]);
mpo_even=cell(1,N);
mpo_odd=cell(1,N);
for j=1:N
    mpo_even{j}=I;
    mpo_odd{j}=I;
end


for j=1:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2;
    LH=LH+1i*H;
    
    L=LH+Lij;
    
    w=expm(dt*L);
    
    if j==1
        w=expm(dt*(L+H1+L1));
    elseif j==(N-1)
        w=expm(dt*(L+HN+LN));
    end
    
    [U,V] = BondMPOsvd(w,dl);
    
    mpo_odd{j}=U;
    mpo_odd{j+1}=V;
    
end
for j=2:2:(N-1)
    
    H=Jx*kron(kron(id,sx),kron(id,sx))+Jy*kron(kron(id,sy),kron(id,sy))...
        +Jz*kron(kron(id,sz),kron(id,sz))...
        +hz(j)*kron(kron(id,sz),Id)/2+hz(j+1)*kron(Id,kron(id,sz))/2;
    LH=-1i*H;
    H=Jx*kron(kron(sx.',id),kron(sx.',id))...
        +Jy*kron(kron(sy.',id),kron(sy.',id))...
        +Jz*kron(kron(sz.',id),kron(sz.',id))...
        +hz(j)*kron(kron(sz,id),Id)/2+hz(j+1)*kron(Id,kron(sz,id))/2;
    LH=LH+1i*H;
    
    L=LH+Lij;
    
    w=expm(dt/2*L);
    
    if j==(N-1)
        w=expm(dt/2*(L+HN+LN));
        
    end
    
    [U,V] = BondMPOsvd(w,dl);
    
    mpo_even{j}=U;
    mpo_even{j+1}=V;
end

MPO=Mpo2Mpo1(mpo_odd,mpo_even,N);
MPO=Mpo2Mpo1(mpo_even,MPO,N);

VMPO=VecMPO(MPO,N);

[MPO,K]=mpsTruncSvd(VMPO,Dmax,EpsilonMPO);
MPO{1}=sqrt(K)*MPO{1};

[VMPO,K]=CanonizeA(VMPO);
VMPO{1}=sqrt(K)*VMPO{1};

[MPO]=VarCompressMPS(VMPO,MPO,precision);
MPO{1}=sqrt(K)*MPO{1};

MPO=MatMPS(MPO,N);