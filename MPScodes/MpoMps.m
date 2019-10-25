function mpsB=MpoMps(mpo,mps,N)
mpsB=cell(1,N);
for i=1:N
     result = contracttensors(mps{i},3,3,mpo{i},4,4);
s_r = size(result);
mpsB{i} = reshape(permute(result,[3 1 4 2 5]),[s_r(1)*s_r(3),s_r(2)*s_r(4),s_r(5)]);
end