
VSCurrent=[];
Nnodes=10;
for node=1:Nnodes

fname = sprintf('MPSJz%dW%dG%dD%dN%dnode%d.mat',Jz,Wd,Gamma,D,N,node);
load(fname)

VSCurrent=[VSCurrent SCurrent'];

hold on
errorbar(1:1:N,Prof,ErrR,'o-')

end

current=mean(real(VSCurrent));
Err=std(real(VSCurrent),0,2)/sqrt(Nnodes*Nrea);

