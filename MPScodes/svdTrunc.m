function [U,S,V]=svdTrunc(T,Dmax,Epsilon)
[m,n]=size(T);
Delt=10^(-13);
if m>=n,
    try
        [U,S,V]=svd(T,0);
    catch
        [U,S,V]=svd(T+Delt*randn(m,n),0);
    end
    Vlam=(diag(S).^2);
Vlam=Vlam/sum(Vlam);
[~,TruncD]=min(abs(cumsum(Vlam)-(1-Epsilon)));
if TruncD > Dmax
    TruncD = Dmax;
end
% Do trimming
U = U(:,1:TruncD);
S = S(1:TruncD,1:TruncD);
S = S/norm(diag(S));
V = V(:,1:TruncD);
else
    try
        [V,S,U]=svd(T',0);
    catch
        [V,S,U]=svd((T+Delt*randn(m,n))',0);
    end
    Vlam=(diag(S).^2);
Vlam=Vlam/sum(Vlam);
[~,TruncD]=min(abs(cumsum(Vlam)-(1-Epsilon)));
if TruncD > Dmax
    TruncD = Dmax;
end
% Do trimming
U = U(1:TruncD,:);
S = S(1:TruncD,1:TruncD);
S = S/norm(diag(S));
V = V(1:TruncD,:);
    
end


V=V';
