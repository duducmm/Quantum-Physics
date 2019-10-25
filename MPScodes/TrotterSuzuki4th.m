function mps=TrotterSuzuki4th(MPO1,MPO2,MPO3,mps,N,D,Epsilon,precision,Sweepmax)


mps=Trotter2ond(MPO1,mps,N,D,Epsilon,precision,Sweepmax);

mps=Trotter2ond(MPO2,mps,N,D,Epsilon,precision,Sweepmax);

mps=Trotter2ond(MPO3,mps,N,D,Epsilon,precision,Sweepmax);

mps=Trotter2ond(MPO2,mps,N,D,Epsilon,precision,Sweepmax);

mps=Trotter2ond(MPO1,mps,N,D,Epsilon,precision,Sweepmax);

