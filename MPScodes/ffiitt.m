function Dist = ffiitt(perfil,II,alpha)

y=alpha(2)*II.^alpha(1);
Dist = sum(abs(perfil-y))*10;
