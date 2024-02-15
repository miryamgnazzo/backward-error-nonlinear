%
% Hadeler problem from the NLEVP collection
%

n = 8;
[F, f] = nlevp('hadeler', n);
[V, L] = be_newton(F, f, [0.217, 0.885]);

