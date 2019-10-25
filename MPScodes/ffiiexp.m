function Dist = ffiiexp(perfil,II,alpha)

y=alpha(1).^(II*alpha(2));
Dist = sum(abs(perfil-y))*10;
