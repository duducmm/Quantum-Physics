function MPO=Mpo2Mpo1(mpo2,mpo1,N)
MPO=cell(1,N);
for i=1:N
     result = contracttensors(mpo1{i},4,3,mpo2{i},4,4);
s_r = size(result);
MPO{i} = reshape(permute(result,[4 1 5 2 6 3]),[s_r(1)*s_r(4),s_r(2)*s_r(5),s_r(6),s_r(3)]);
end