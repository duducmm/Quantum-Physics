function VMPO=VecMPO(MPO,N)
VMPO=cell(1,N);
for i=1:N
s_r = size(MPO{i});
VMPO{i} = reshape(MPO{i},[s_r(1),s_r(2),s_r(3)*s_r(4)]);
end