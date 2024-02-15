%
% Generate a pencil
%

n = 12;
k = 2; % Number of f_j(z)
p = 3; % Number of eigenvalues to consider

epsilon = 1e-3;

A = randn(n, n); A = A + A';
B = randn(n, n); B = B * B';

dA = randn(n, n); dA = dA + dA'; dA = dA / norm(dA, 'fro') * epsilon;
dB = randn(n, n); dB = dB * dB'; dB = dB / norm(dB, 'fro') * epsilon;

Ap = A + dA;
Bp = B + dB;

[V, L] = eig(Ap, -Bp);
V = V(:, 1:p); L = L(1:p, 1:p);

% norm(Ap * V + Bp * V * L, 'fro')

Res = A * V + B * V * L;
norm(Res, 'fro')

[Q, R] = qr(V);
Qt  = Q(:, 1:p);
Rt = R(1:p, :);
Qto = Q(:, p+1:end);

A21 = Qto' * A * Qt;
B21 = Qto' * B * Qt;

Res1 = - Qt' * Res;
Res2 = - Qto' * Res;

R1 = Rt;
R2 = Rt * L;

% Solve for A21, B21
Block_21 = Res2 * pinv([ R1 ; R2 ]);
A21t = Block_21(:, 1:p);
B21t = Block_21(:, p+1:2*p);

% Solve for the symmetric blocks A11, B11
P = 1 : p^2; P = reshape(P, p, p); P = P'; P = P(:);
CM = eye(p^2); CM = CM(:, P);

Sys_11 = [ ...
    kron([ R1', R2' ], eye(p)) ; ...
    blkdiag(CM - eye(p^2), CM - eye(p^2))
];

Block_11 = pinv(Sys_11) * [ Res1(:) ; zeros(2*p^2, 1) ];

A11t = reshape(Block_11(1:p^2), p, p);
B11t = reshape(Block_11(p^2+1:end), p, p);

At = [ A11t, A21t' ; A21t , zeros(n-p, n-p) ];
Bt = [ B11t, B21t' ; B21t , zeros(n-p, n-p) ];

Ab = Q * At * Q';
Bb = Q * Bt * Q';

Resb = (A + Ab) * V + (B + Bb) * V * L;
norm(Resb, 'fro')
