%Embedding existing skyrmion into wider lattice
%Needs loaded configuration
%Advisable to clear before
newsize = 71;
oldsize = 41;
delta = abs(newsize - oldsize)/2;
left_b = (delta)+1;
right_b = newsize-delta;

Tx = zeros(newsize,newsize);
Ty = Tx;
Tz = ones(newsize,newsize);
Tx(left_b:right_b,left_b:right_b) = Sx;
Ty(left_b:right_b,left_b:right_b) = Sy;
Tz(left_b:right_b,left_b:right_b) = Sz;
Sx = Tx;
Sy = Ty;
Sz = Tz;
