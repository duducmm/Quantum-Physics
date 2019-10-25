function MPO=MatMPS(MPS,N)
MPO=cell(1,N);
for i=1:N
s_r = size(MPS{i});
MPO{i} = reshape(MPS{i},[s_r(1),s_r(2),sqrt(s_r(3)),sqrt(s_r(3))]);
%MPO{i}=permute(MPO{i},[1 2 4 3]);
end