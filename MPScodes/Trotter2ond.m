function mps=Trotter2ond(MPO,mps,N,D,Epsilon,precision,Sweepmax)

mps=MpoMps(MPO,mps,N);
%mps=CanonizeA(mps);
[mpsSVD]=mpsTruncSvd2(mps,D,Epsilon);
mpsSVD=CanonizeA(mpsSVD);
[mps]=VarCompressMPS(mps,mpsSVD,precision,Sweepmax);